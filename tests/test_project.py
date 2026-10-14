"""Synthetic test fixtures verify implementation; they are not research data."""
import json
from pathlib import Path
import sys
import subprocess

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core import (contact_features, fixed_rmsd, grouped_predictions,
                  paired_bootstrap, rank_metrics, validate_xyz)
from pipeline import analyze, download, figures, pdb_atoms, report
from dataset import validate_structure, validate_chemistry, prepare_data


def test_rmsd_preserves_placement_error():
    a = np.array([[0., 0., 0.], [1., 0., 0.]])
    assert fixed_rmsd(a, a + [3, 0, 0]) == 3
    assert fixed_rmsd(a, a) == 0


@pytest.mark.parametrize('x', [[], [[1, 2]], [[float('nan'), 0, 0]], [[float('inf'), 0, 0]]])
def test_invalid_coordinates(x):
    with pytest.raises(ValueError):
        validate_xyz(x)


def test_rmsd_atom_count_mismatch():
    with pytest.raises(ValueError):
        fixed_rmsd([[0, 0, 0]], [[0, 0, 0], [1, 0, 0]])


def test_contacts_rigid_transform_invariant():
    l = np.array([[0., 0, 0], [1, 2, 0]])
    p = np.array([[3., 0, 0], [2, 3, 1]])
    rotation = np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]])
    a = contact_features(l, ['C', 'N'], p, ['C', 'O'], -4)
    b = contact_features(l @ rotation + 10, ['C', 'N'], p @ rotation + 10, ['C', 'O'], -4)
    np.testing.assert_allclose(a, b)


def test_clash_descriptor_penalizes_overlap():
    close = contact_features([[0, 0, 0]], ['C'], [[1, 0, 0]], ['C'], -3)
    far = contact_features([[0, 0, 0]], ['C'], [[5, 0, 0]], ['C'], -3)
    assert close[1] > far[1]


def test_contact_metadata_validation():
    with pytest.raises(ValueError):
        contact_features([[0, 0, 0]], [], [[1, 0, 0]], ['C'], -3)


def fixture():
    rng = np.random.default_rng(41)
    x = rng.normal(size=(40, 7))
    y = np.abs(rng.normal(size=40)) * 3
    groups = np.array(['s1:a'] * 10 + ['s1:b'] * 10 + ['s2:c'] * 10 + ['s2:d'] * 10)
    return x, y, groups


def test_no_held_series_label_leakage():
    x, y, g = fixture()
    a, folds = grouped_predictions(x, y, g)
    changed = y.copy()
    changed[g == 's2:c'] += 100
    b, _ = grouped_predictions(x, changed, g)
    np.testing.assert_array_equal(a[np.char.startswith(g, 's2')], b[np.char.startswith(g, 's2')])
    for fold in folds:
        assert set(fold['train_fragments']).isdisjoint(fold['test_fragments'])
        assert all(z.startswith(fold['held_out_series']) for z in fold['test_fragments'])


def test_shuffled_control_reproducible():
    x, y, g = fixture()
    a, _ = grouped_predictions(x, y, g, shuffle=True)
    b, _ = grouped_predictions(x, y, g, shuffle=True)
    np.testing.assert_array_equal(a, b)


def test_one_series_rejected():
    with pytest.raises(ValueError):
        grouped_predictions(np.ones((4, 7)), np.ones(4), ['s1:a'] * 4)


def test_ranking_and_oracle_distinguish_sampling_failure():
    a = rank_metrics([4, 1, 3], [-5, -4, -3])
    assert a['top1_success'] == 0
    assert a['top3_success'] == 1
    assert a['oracle_success'] == 1
    b = rank_metrics([4, 5, 3], [-5, -4, -3])
    assert b['oracle_success'] == 0


def test_threshold_and_ties():
    assert rank_metrics([2, 5], [0, 0])['top1_success'] == 1
    assert rank_metrics([5, 2], [0, 0])['top1_success'] == 0


def test_paired_bootstrap_constant_difference():
    result = paired_bootstrap([1, 2, 3], [2, 3, 4])
    assert result['mean_difference'] == -1
    assert result['percentile_95_CI'] == [-1, -1]
    assert result['n_fragments'] == 3


def test_offline_missing_data_fails(tmp_path):
    with pytest.raises(FileNotFoundError):
        download('https://example.invalid/data', tmp_path / 'missing.pdb', True)


def test_invalid_cache_rejected(tmp_path):
    p = tmp_path / 'bad.pdb'
    p.write_text('oops')
    with pytest.raises(ValueError):
        download('https://example.invalid/data', p, True)


def test_parser_occupancy_altloc_and_first_model(tmp_path):
    def line(serial, alt, x, occ):
        return f'ATOM  {serial:5d}  CA {alt}ALA A   1    {x:8.3f}{0:8.3f}{0:8.3f}{occ:6.2f}{20:6.2f}          C '
    p = tmp_path / 'fixture.pdb'
    p.write_text('\n'.join([line(1, 'A', 0, .4), line(2, 'B', 1, .6), 'ENDMDL', line(3, ' ', 9, 1)]))
    atoms = pdb_atoms(p)
    assert len(atoms) == 1
    assert atoms[0]['xyz'][0] == 1
    assert atoms[0]['line'][16] == ' '


def test_evaluation_artifacts(tmp_path):
    x, y, g = fixture()
    records = [{'pdb': name.split(':')[1], 'pose': i} for i, name in enumerate(g)]
    metrics, predictions, summary = analyze(x, y, g, records, tmp_path, 17)
    report(metrics, summary, tmp_path)
    assert np.isfinite(predictions).all()
    assert len(metrics) == 12
    assert json.loads((tmp_path / 'summary.json').read_text())['independent_series'] == 2
    assert 'no MD' in (tmp_path / 'REPORT.md').read_text()
    assert (tmp_path / 'poses.csv').read_text().count('\n') == 41


def test_figure_generation(tmp_path):
    class Conformer:
        def GetPositions(self):
            return np.array([[0., 0, 0], [1, 0, 0], [1, 1, 0]])
    class Bond:
        def GetBeginAtomIdx(self):
            return 0
        def GetEndAtomIdx(self):
            return 1
    class Molecule:
        def GetConformer(self):
            return Conformer()
        def GetBonds(self):
            return [Bond()]
    x, y, g = fixture()
    rows = [{'pdb': name.split(':')[1], 'pose': i} for i, name in enumerate(g)]
    metrics, _, _ = analyze(x, y, g, rows, tmp_path, 17)
    protein = [{'name': 'CA', 'xyz': [2, 3, 4]}, {'name': 'CA', 'xyz': [3, 4, 5]}]
    figures(x, y, g, metrics, [('TEST FIXTURE', Molecule(), protein, Molecule())], tmp_path)
    for name in ('ranking.png', 'poses.png'):
        assert (tmp_path / name).read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    assert '<svg' in (tmp_path / 'ranking.svg').read_text()


def test_dependency_manifest_complete_and_pinned():
    root = Path(__file__).resolve().parents[1]
    requirements = [l for l in (root / 'requirements.txt').read_text().splitlines()
                    if l.strip() and not l.startswith('#')]
    names = {l.split('==')[0] for l in requirements}
    assert names == {'numpy', 'scipy', 'scikit-learn', 'matplotlib', 'rdkit',
                     'biopython', 'meeko', 'gemmi', 'pytest'}
    assert all(len(l.split('==')) == 2 for l in requirements)
    assert not (root / 'environment.yml').exists()


def test_launcher_shell_syntax():
    script = Path(__file__).resolve().parents[1] / 'run.sh'
    subprocess.run(['bash', '-n', str(script)], check=True)


def test_dataset_rejects_wrong_structure_identity(tmp_path):
    p = tmp_path / 'wrong.pdb'
    p.write_text('HEADER'.ljust(62) + 'XXXX\nEXPDTA    X-RAY DIFFRACTION\n')
    with pytest.raises(ValueError, match='identity'):
        validate_structure(p, '8G65', 'YW9')


def test_dataset_rejects_missing_coordinates(tmp_path):
    p = tmp_path / 'empty.pdb'
    p.write_text('HEADER'.ljust(62) + '8G65\nEXPDTA    X-RAY DIFFRACTION\n')
    with pytest.raises(ValueError, match='truncated'):
        validate_structure(p, '8G65', 'YW9')


def test_dataset_rejects_wrong_component_identity(tmp_path):
    p = tmp_path / 'wrong.cif'
    p.write_text('data_OTHER\n')
    with pytest.raises(ValueError, match='identity'):
        validate_chemistry(p, 'YW9')


def test_dataset_rejects_incomplete_chemistry(tmp_path):
    p = tmp_path / 'empty.cif'
    p.write_text('data_YW9\n')
    with pytest.raises(ValueError, match='missing chemistry'):
        validate_chemistry(p, 'YW9')


def test_missing_real_data_cannot_produce_dataset_manifest(tmp_path):
    with pytest.raises(FileNotFoundError):
        prepare_data(tmp_path, offline=True)
    assert not (tmp_path / 'manifest.json').exists()


def test_cached_dataset_hash_change_rejected(tmp_path):
    (tmp_path / '8G65.pdb').write_text('fixture bytes only' * 20)
    manifest = {'sources': [{'file': '8G65.pdb', 'sha256': '0' * 64}]}
    (tmp_path / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='changed since'):
        prepare_data(tmp_path, offline=True)
    assert json.loads((tmp_path / 'manifest.json').read_text()) == manifest


def real_files():
    from dataset import ROOT, ENTRIES
    directory = ROOT / 'data'
    if not all((directory / f'{p}.pdb').exists() and (directory / f'{l}.cif').exists()
               for p, l, _ in ENTRIES):
        pytest.skip('Real-data integration: first run python3 dataset.py to download RCSB files')
    return directory, ENTRIES


def test_real_dataset_chemical_graphs_and_chirality():
    directory, entries = real_files()
    from pipeline import structure
    from rdkit import Chem
    for pdb, ligand, _ in entries:
        reference, protein, metadata = structure(directory / f'{pdb}.pdb', directory / f'{ligand}.cif', ligand)
        assert reference.GetNumAtoms() == metadata['heavy_atoms']
        assert len(protein) > 500
        assert np.isfinite(reference.GetConformer().GetPositions()).all()
        assert all(tag != '?' for _, tag in Chem.FindMolChiralCenters(reference, includeUnassigned=True))


def test_real_ligands_meeko_roundtrip():
    directory, entries = real_files()
    from pipeline import structure
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
    for pdb, ligand, _ in entries:
        reference, _, _ = structure(directory / f'{pdb}.pdb', directory / f'{ligand}.cif', ligand)
        start = Chem.AddHs(Chem.Mol(reference))
        start.RemoveAllConformers()
        assert AllChem.EmbedMolecule(start, randomSeed=17) == 0
        setups = MoleculePreparation().prepare(start)
        text, ok, error = PDBQTWriterLegacy.write_string(setups[0])
        assert ok, error
        recovered = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(text, skip_typing=True))
        assert len(recovered) == 1 and recovered[0] is not None
        heavy = Chem.RemoveHs(recovered[0])
        assert heavy.GetNumAtoms() == reference.GetNumAtoms()
        assert heavy.HasSubstructMatch(reference, useChirality=True)


def test_real_vina_smoke_run(tmp_path):
    directory, entries = real_files()
    from pipeline import dock, structure
    from rdkit.Chem import rdMolAlign
    pdb, ligand, _ = entries[0]
    ref, protein, _ = structure(directory / f'{pdb}.pdb', directory / f'{ligand}.cif', ligand)
    poses, scores = dock(ref, protein, tmp_path, [17], exhaustiveness=1, cpu=1)
    assert len(poses) == len(scores) and len(poses) > 0
    assert np.isfinite(scores).all()
    assert all(np.isfinite(rdMolAlign.CalcRMS(p, ref)) for p in poses)
