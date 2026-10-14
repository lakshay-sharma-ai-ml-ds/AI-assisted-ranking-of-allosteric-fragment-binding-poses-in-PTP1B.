"""Experimental PTP1B redocking and out-of-series pose rescoring.

This is a retrospective, known-pocket study. No native pose is inserted into
the candidate pool. Docking starts from independently embedded coordinates.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys

import numpy as np
from core import FEATURES, contact_features, grouped_predictions, paired_bootstrap, rank_metrics

from dataset import ENTRIES, PAPER, download, prepare_data


def pdb_atoms(path):
    """First model, highest-occupancy alternate per atom; heavy atoms only."""
    atoms = {}
    for line in path.read_text().splitlines():
        if line.startswith('ENDMDL'):
            break
        if line[:6].strip() not in ('ATOM', 'HETATM'):
            continue
        element = line[76:78].strip().upper() or line[12:16].strip()[0]
        if element in ('H', 'D'):
            continue
        key = (line[:6].strip(), line[21], line[22:27], line[17:20].strip(), line[12:16].strip())
        a = {'record': key[0], 'chain': key[1], 'resid': key[2], 'resname': key[3],
             'name': key[4], 'xyz': [float(line[30:38]), float(line[38:46]), float(line[46:54])],
             'occupancy': float(line[54:60].strip() or 0), 'element': element,
             'bfactor': float(line[60:66].strip() or 0), 'line': line[:16] + ' ' + line[17:]}
        if key not in atoms or a['occupancy'] > atoms[key]['occupancy']:
            atoms[key] = a
    if not atoms:
        raise ValueError(f'No atoms in {path}')
    return list(atoms.values())


def structure(path, ccd, ligand_id):
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from rdkit import Chem
    from scipy.spatial.distance import cdist
    d = MMCIF2Dict(str(ccd))
    names = d['_chem_comp_atom.atom_id']
    symbols = d['_chem_comp_atom.type_symbol']
    heavy = [n for n, e in zip(names, symbols) if e not in ('H', 'D')]
    atoms = pdb_atoms(path)
    ligand_residues = sorted({(a['chain'], a['resid']) for a in atoms if a['resname'] == ligand_id})
    candidates = []
    chains = sorted({a['chain'] for a in atoms if a['record'] == 'ATOM'})
    for chain, resid in ligand_residues:
        lig = {a['name']: a for a in atoms if a['resname'] == ligand_id and
               a['chain'] == chain and a['resid'] == resid}
        if not all(n in lig for n in heavy):
            continue
        for pc in chains:
            protein = [a for a in atoms if a['record'] == 'ATOM' and a['chain'] == pc]
            dist = cdist(np.array([lig[n]['xyz'] for n in heavy]), np.array([a['xyz'] for a in protein]))
            contact = int(np.sum(dist < 4.5))
            candidates.append((contact, chain, resid, pc, lig, protein))
    if not candidates:
        raise ValueError(f'No complete {ligand_id} instance in {path}; check current CCD identity')
    # Deterministic contact-based choice; other chains are excluded from monomer receptor.
    _, lc, resid, pc, lig, protein = max(candidates, key=lambda a: (a[0], a[1], a[2], a[3]))
    rw = Chem.RWMol()
    index = {}
    charges = d['_chem_comp_atom.charge']
    for name, symbol, charge in zip(names, symbols, charges):
        if name not in heavy:
            continue
        atom = Chem.Atom(symbol.capitalize())
        atom.SetFormalCharge(int(charge) if charge not in ('.', '?') else 0)
        index[name] = rw.AddAtom(atom)
    orders = {'SING': Chem.BondType.SINGLE, 'DOUB': Chem.BondType.DOUBLE,
              'TRIP': Chem.BondType.TRIPLE, 'AROM': Chem.BondType.AROMATIC}
    for a, b, order in zip(d['_chem_comp_bond.atom_id_1'], d['_chem_comp_bond.atom_id_2'],
                            d['_chem_comp_bond.value_order']):
        if a in index and b in index:
            rw.AddBond(index[a], index[b], orders[order])
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    conf = Chem.Conformer(len(index))
    for name, i in index.items():
        conf.SetAtomPosition(i, lig[name]['xyz'])
    mol.AddConformer(conf)
    Chem.AssignAtomChiralTagsFromStructure(mol, confId=0, replaceExistingTags=True)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    stereo = dict(zip(names, d.get('_chem_comp_atom.pdbx_stereo_config', ['N'] * len(names))))
    for name, idx in index.items():
        if stereo[name] in ('R', 'S'):
            atom = mol.GetAtomWithIdx(idx)
            if not atom.HasProp('_CIPCode') or atom.GetProp('_CIPCode') != stereo[name]:
                raise ValueError(f'Crystal/CCD stereochemistry mismatch at {name}')
    if any(tag == '?' for _, tag in Chem.FindMolChiralCenters(mol, includeUnassigned=True)):
        raise ValueError('Unassigned ligand stereochemistry')
    metadata = {'ligand': ligand_id, 'ligand_chain': lc, 'ligand_residue': resid.strip(),
                'protein_chain': pc, 'heavy_atoms': len(index),
                'mean_occupancy': float(np.mean([lig[n]['occupancy'] for n in heavy])),
                'mean_B_A2': float(np.mean([lig[n]['bfactor'] for n in heavy]))}
    return mol, protein, metadata


def command(args, logfile):
    result = subprocess.run([str(a) for a in args], capture_output=True, text=True)
    logfile.write_text(result.stdout + '\n' + result.stderr)
    if result.returncode:
        raise RuntimeError(f'{args[0]} failed; inspect {logfile}')


def dock(reference, protein, folder, seeds, exhaustiveness, cpu):
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
    from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
    folder.mkdir(parents=True, exist_ok=True)
    protein_pdb = folder / 'receptor.pdb'
    protein_pdb.write_text('\n'.join(a['line'] for a in protein) + '\nTER\nEND\n')
    receptor = folder / 'receptor.pdbqt'
    command(['obabel', '-ipdb', protein_pdb, '-opdbqt', '-O', receptor,
             '-xr', '-h', '--partialcharge', 'gasteiger'], folder / 'receptor-preparation.log')
    text = receptor.read_text()
    if 'ATOM' not in text or 'BRANCH' in text:
        raise ValueError('Invalid rigid receptor PDBQT')
    receptor_heavy = [l for l in text.splitlines() if l[:6].strip() in ('ATOM', 'HETATM')
                      and l.split()[-1] not in ('H', 'HD', 'HS')]
    if len(receptor_heavy) != len(protein):
        raise ValueError('Receptor preparation changed the heavy-atom count; inspect preparation log')
    center = np.mean(reference.GetConformer().GetPositions(), axis=0)
    extent = np.ptp(reference.GetConformer().GetPositions(), axis=0)
    box = np.maximum(extent + 12, 20)
    settings = {'center_A': center.tolist(), 'box_size_A': box.tolist(), 'seeds': seeds,
                'exhaustiveness': exhaustiveness, 'cpu': cpu, 'candidate_source': 'AutoDock Vina',
                'native_pose_injected': False, 'receptor': 'native holo monomer, rigid'}
    poses, scores = [], []
    for seed in seeds:
        # Crystal coordinates are discarded; ETKDG starts from chemistry alone.
        start = Chem.AddHs(Chem.Mol(reference))
        start.RemoveAllConformers()
        p = AllChem.ETKDGv3()
        p.randomSeed = seed
        if AllChem.EmbedMolecule(start, p) != 0:
            raise ValueError('Independent ligand embedding failed')
        if not AllChem.MMFFHasAllMoleculeParams(start):
            raise ValueError('Missing MMFF ligand parameters')
        AllChem.MMFFOptimizeMolecule(start)
        setup = MoleculePreparation().prepare(start)
        if len(setup) != 1:
            raise ValueError('Unexpected multiple ligand setups')
        ligand_text, ok, error = PDBQTWriterLegacy.write_string(setup[0])
        if not ok:
            raise ValueError(error)
        ligand = folder / f'ligand-{seed}.pdbqt'
        ligand.write_text(ligand_text)
        output = folder / f'poses-{seed}.pdbqt'
        args = ['vina', '--receptor', receptor, '--ligand', ligand, '--out', output,
                '--seed', seed, '--cpu', cpu, '--exhaustiveness', exhaustiveness,
                '--num_modes', 20, '--energy_range', 10]
        for axis, c, s in zip('xyz', center, box):
            args += [f'--center_{axis}', float(c), f'--size_{axis}', float(s)]
        command(args, folder / f'vina-{seed}.log')
        output_text = output.read_text()
        energy = [float(l.split()[3]) for l in output_text.splitlines() if l.startswith('REMARK VINA RESULT:')]
        exported = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(output_text, skip_typing=True))
        if len(exported) != 1 or exported[0] is None:
            raise ValueError('Could not reconstruct docked ligand with atom identities')
        heavy = Chem.RemoveHs(exported[0])
        if heavy.GetNumConformers() != len(energy):
            raise ValueError('Pose/energy count mismatch')
        if not heavy.HasSubstructMatch(Chem.RemoveHs(reference), useChirality=True) or heavy.GetNumAtoms() != reference.GetNumAtoms():
            raise ValueError('Docked/reference chemical graphs differ')
        for conf, score in zip(heavy.GetConformers(), energy):
            single = Chem.Mol(heavy)
            single.RemoveAllConformers()
            single.AddConformer(Chem.Conformer(conf), assignId=True)
            poses.append(single)
            scores.append(score)
    # Merge replicate searches, retaining the lower-scoring duplicate.
    selected, selected_scores = [], []
    for i in np.argsort(scores):
        if all(rdMolAlign.CalcRMS(poses[i], q) >= 0.5 for q in selected):
            selected.append(poses[i])
            selected_scores.append(scores[i])
    if not selected:
        raise ValueError('Vina returned no valid candidate poses')
    with Chem.SDWriter(str(folder / 'candidates.sdf')) as writer:
        for pose, score in zip(selected, selected_scores):
            pose.SetDoubleProp('vina_score_kcal_mol', float(score))
            writer.write(pose)
    (folder / 'settings.json').write_text(json.dumps(settings, indent=2))
    return selected, selected_scores


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze(x, rmsd, groups, records, out, seed):
    predictions, folds = grouped_predictions(x, rmsd, groups, seed)
    shuffled, _ = grouped_predictions(x, rmsd, groups, seed, shuffle=True)
    metrics = []
    for g in sorted(set(groups)):
        sel = np.array(groups) == g
        for method, score in [('Vina', x[:, 0]), ('ML', predictions), ('Shuffled-label ML', shuffled)]:
            metrics.append({'fragment': g.split(':')[1], 'series': g.split(':')[0], 'method': method,
                            'n_poses': int(sel.sum()), **rank_metrics(rmsd[sel], score[sel])})
    for row, pred, control in zip(records, predictions, shuffled):
        row['out_of_series_predicted_log1p_rmsd'] = float(pred)
        row['shuffled_label_prediction'] = float(control)
    write_csv(out / 'poses.csv', records)
    write_csv(out / 'metrics.csv', metrics)
    a = [m['top1_rmsd_A'] for m in metrics if m['method'] == 'ML']
    b = [m['top1_rmsd_A'] for m in metrics if m['method'] == 'Vina']
    summary = {'paired_ML_minus_Vina_top1_rmsd_A': paired_bootstrap(a, b, seed), 'folds': folds,
               'success_threshold_A': 2, 'independent_series': 2,
               'warning': 'Five related fragments in two series; bootstrap CI is descriptive, not confirmatory.'}
    (out / 'summary.json').write_text(json.dumps(summary, indent=2))
    return metrics, predictions, summary


def figures(x, rmsd, groups, metrics, examples, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'savefig.dpi': 180})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ids = sorted({m['fragment'] for m in metrics})
    for j, (method, color) in enumerate([('Vina', '#64748b'), ('ML', '#0f766e')]):
        vals = [next(m['top1_rmsd_A'] for m in metrics if m['fragment'] == p and m['method'] == method) for p in ids]
        axes[0].bar(np.arange(len(ids)) + j * .35, vals, .35, label=method, color=color)
    axes[0].axhline(2, ls='--', color='#b45309', label='2 Å success threshold')
    axes[0].set(xticks=np.arange(len(ids)) + .175, xticklabels=ids, ylabel='Top-ranked symmetry RMSD (Å)', title='Held-out chemical-series ranking')
    axes[0].legend(fontsize=8)
    for g in sorted(set(groups)):
        sel = np.array(groups) == g
        axes[1].scatter(x[sel, 0], rmsd[sel], s=16, alpha=.7, label=g.split(':')[1])
    axes[1].set(xlabel='Vina score (kcal/mol)', ylabel='Crystal-pose RMSD (Å)', title='Docking score versus pose correctness')
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / 'ranking.png')
    fig.savefig(out / 'ranking.svg')
    plt.close(fig)
    fig = plt.figure(figsize=(13, 8))
    for i, (pdb, reference, protein, pose) in enumerate(examples):
        ax = fig.add_subplot(2, 3, i + 1, projection='3d')
        ca = np.array([a['xyz'] for a in protein if a['name'] == 'CA'])
        native = reference.GetConformer().GetPositions()
        predicted = pose.GetConformer().GetPositions()
        center = native.mean(axis=0)
        ca = ca[np.linalg.norm(ca - center, axis=1) < 14]
        ax.scatter(*ca.T, s=5, color='#cbd5e1', alpha=.5)
        ax.scatter(*native.T, s=30, color='#0f766e', label='Crystal')
        ax.scatter(*predicted.T, s=23, color='#d97706', label='ML top-ranked')
        for mol, color in [(reference, '#0f766e'), (pose, '#d97706')]:
            xyz = mol.GetConformer().GetPositions()
            for bond in mol.GetBonds():
                a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                ax.plot(*xyz[[a, b]].T, color=color, lw=1.5)
        for axis in 'xyz':
            getattr(ax, f'set_{axis}lim')(center['xyz'.index(axis)] - 10, center['xyz'.index(axis)] + 10)
        ax.set_box_aspect((1, 1, 1))
        ax.set_title(pdb)
        ax.set(xlabel='x (Å)', ylabel='y (Å)', zlabel='z (Å)')
    ax.legend(fontsize=8)
    fig.suptitle('Reference and predicted poses in the receptor coordinate frame')
    fig.tight_layout()
    fig.savefig(out / 'poses.png')
    plt.close(fig)


def report(metrics, summary, out):
    lines = ['# PTP1B pose-ranking results', '', 'Retrospective known-pocket holo redocking; no MD performed.', '',
             '| Fragment | Method | Top-1 RMSD (Å) | Top-1 success | Top-3 success | Oracle RMSD (Å) |',
             '|---|---|---:|---:|---:|---:|']
    for m in metrics:
        lines.append(f"| {m['fragment']} | {m['method']} | {m['top1_rmsd_A']:.3f} | {m['top1_success']} | {m['top3_success']} | {m['oracle_rmsd_A']:.3f} |")
    delta = summary['paired_ML_minus_Vina_top1_rmsd_A']
    lines += ['', f"Mean paired ML − Vina top-1 RMSD: {delta['mean_difference']:.3f} Å. Negative favors ML.",
              f"Descriptive fragment-bootstrap 95% interval: {delta['percentile_95_CI']}. Two independent chemical series only.",
              '', '![Ranking](ranking.png)', '', '![Pose overlays](poses.png)', '',
              '## Theory and interpretation', '',
              'RMSD is the square root of the mean squared heavy-atom displacement, minimized over graph automorphisms.',
              'Ligands are never independently aligned to the crystal ligand: that would erase placement errors.',
              'Vina samples conformations inside a known pocket and scores approximate interactions. Its scores are not measured binding affinities.',
              'A standardized ridge model (alpha=10) predicts log(1+RMSD) from Vina score and receptor-contact features.',
              'No native coordinates, native contacts, atom displacement, pocket ID, or fragment identity enter its feature vector.',
              'The entire held-out chemical series is excluded from training. The two folds test transfer between different allosteric pockets as well as chemistry.',
              'The candidate oracle separates sampling failure from ranking failure; shuffled training labels provide a negative control.',
              'Polar proximity is a coarse contact descriptor, not an angle-aware hydrogen-bond calculation.',
              '', '## Limits', '',
              'Five structures are insufficient to establish generalizable AI performance. Poses are not independent experiments.',
              'Holo receptors and crystal-derived box centers make this an optimistic retrospective benchmark; it does not discover cryptic pockets.',
              'One selected monomer is used. Alternate conformations, crystal packing, missing residues, protonation, waters, and receptor flexibility can affect results.',
              'Inspect data-quality.json and the raw preparation logs before interpreting a result. Default ligand protonation follows the deposited CCD graph.',
              'This workflow does not reproduce the Anton study, run MD, establish affinity, or identify a new drug.',
              '', f'Paper: {PAPER}']
    (out / 'REPORT.md').write_text('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true', help='Require cached structural files')
    parser.add_argument('--download-only', action='store_true')
    parser.add_argument('--seeds', nargs='+', type=int, default=[17, 41, 83])
    parser.add_argument('--exhaustiveness', type=int, default=32)
    parser.add_argument('--cpu', type=int, default=1)
    args = parser.parse_args()
    if args.cpu < 1 or args.exhaustiveness < 1 or any(s <= 0 or s >= 2**31 for s in args.seeds):
        parser.error('CPU, exhaustiveness and seeds must be positive; seeds must fit a signed 32-bit integer')
    root = Path(__file__).resolve().parent
    data, out = root / 'data', root / 'results'
    out.mkdir(exist_ok=True)
    # Never leave a prior successful report looking like this run's output.
    if not args.download_only:
        for filename in ('REPORT.md', 'summary.json', 'manifest.json'):
            (out / filename).unlink(missing_ok=True)
    sources = prepare_data(data, args.offline)
    if args.download_only:
        print('Downloaded five experimental structures and chemical-component definitions.')
        return
    for binary in ('obabel', 'vina'):
        if not shutil.which(binary):
            raise RuntimeError(f'Missing {binary}; use bash run.sh to create the conda environment')
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    x, y, groups, records, objects, quality = [], [], [], [], [], []
    for pdb, ligand, series in ENTRIES:
        print(f'Preparing and docking {pdb} ({ligand})', flush=True)
        reference, protein, metadata = structure(data / f'{pdb}.pdb', data / f'{ligand}.cif', ligand)
        quality.append({'pdb': pdb, **metadata})
        poses, scores = dock(reference, protein, out / pdb, args.seeds, args.exhaustiveness, args.cpu)
        objects.append((pdb, reference, protein, poses))
        protein_xyz = np.array([a['xyz'] for a in protein])
        protein_elements = [a['element'] for a in protein]
        for i, (pose, score) in enumerate(zip(poses, scores)):
            features = contact_features(pose.GetConformer().GetPositions(),
                                        [a.GetSymbol().upper() for a in pose.GetAtoms()],
                                        protein_xyz, protein_elements, score)
            rmsd = float(rdMolAlign.CalcRMS(pose, reference))
            x.append(features)
            y.append(rmsd)
            groups.append(f'{series}:{pdb}')
            records.append({'pdb': pdb, 'ligand': ligand, 'series': series, 'pose': i,
                            'crystal_rmsd_A': rmsd, **dict(zip(FEATURES, map(float, features)))})
    x, y = np.array(x), np.array(y)
    metrics, predictions, summary = analyze(x, y, groups, records, out, args.seeds[0])
    examples = []
    for pdb, ref, protein, poses in objects:
        idx = [i for i, r in enumerate(records) if r['pdb'] == pdb]
        top = int(np.argmin(predictions[idx]))
        examples.append((pdb, ref, protein, poses[top]))
        with Chem.SDWriter(str(out / pdb / 'ML-top.sdf')) as writer:
            writer.write(poses[top])
        with Chem.SDWriter(str(out / pdb / 'crystal-reference.sdf')) as writer:
            writer.write(ref)
    figures(x, y, groups, metrics, examples, out)
    report(metrics, summary, out)
    versions = {}
    from importlib.metadata import version
    for pkg in ['numpy', 'scipy', 'scikit-learn', 'matplotlib', 'rdkit', 'meeko', 'biopython']:
        versions[pkg] = version(pkg)
    manifest = {'paper': PAPER, 'sources': sources, 'parameters': vars(args),
                'python': sys.version, 'platform': platform.platform(), 'versions': versions,
                'code_sha256': {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                                for name in ('pipeline.py', 'core.py', 'dataset.py', 'datasets.csv')},
                'protocol': 'known-pocket holo redocking, leave-chemical-series-out ML; no MD'}
    for binary in ('vina', 'obabel'):
        v = subprocess.run([binary, '--version' if binary == 'vina' else '-V'], capture_output=True, text=True)
        manifest[binary] = v.stdout + v.stderr
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    (out / 'data-quality.json').write_text(json.dumps(quality, indent=2))
    print(f'Finished. Open {out / "REPORT.md"}')


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(1)
