"""Download and validate the experimental dataset using Python's standard library.

No generated structures, artificial fallback data, or package installation.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time
import urllib.request

ROOT = Path(__file__).resolve().parent
PAPER = 'https://doi.org/10.1021/acs.jcim.3c00236'
with (ROOT / 'datasets.csv').open(newline='') as handle:
    DATASETS = list(csv.DictReader(handle))
ENTRIES = [(r['pdb'], r['ligand'], r['series']) for r in DATASETS]


def download(url, path, offline):
    if path.exists():
        if path.stat().st_size < 100:
            raise ValueError(f'Invalid cached file: {path}')
        return
    if offline:
        raise FileNotFoundError(f'Missing {path}; download {url}')
    path.parent.mkdir(parents=True, exist_ok=True)
    last = None
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'PTP1B-research/1.0'})
            with urllib.request.urlopen(request, timeout=30) as response:
                data = response.read()
            if len(data) < 100 or data.lstrip().lower().startswith((b'<!doctype html', b'<html')):
                raise ValueError('Server returned invalid structural data')
            temp = path.with_suffix(path.suffix + '.tmp')
            temp.write_bytes(data)
            temp.replace(path)
            return
        except (OSError, ValueError) as exc:
            last = exc
            if attempt < 2:
                time.sleep(attempt + 1)
    raise RuntimeError(f'Cannot download {url}. Put the file at {path} and use --offline. {last}')


def validate_structure(path, pdb, ligand):
    text = path.read_text()
    lines = text.splitlines()
    headers = [l for l in lines if l.startswith('HEADER')]
    if not headers or headers[0][62:66].strip().upper() != pdb:
        raise ValueError(f'{path.name}: PDB identity does not match {pdb}')
    if not any(l.startswith('EXPDTA') and 'X-RAY DIFFRACTION' in l for l in lines):
        raise ValueError(f'{path.name}: expected an experimental X-ray structure')
    atom_lines = [l for l in lines if l.startswith('ATOM  ')]
    if len(atom_lines) < 500:
        raise ValueError(f'{path.name}: missing/truncated protein coordinates')
    ligand_lines = [l for l in lines if l.startswith('HETATM') and l[17:20].strip() == ligand]
    if len(ligand_lines) < 8:
        raise ValueError(f'{path.name}: expected ligand {ligand} is missing/truncated')
    for l in atom_lines + ligand_lines:
        xyz = [float(l[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))]
        if any(not (-100000 < x < 100000) for x in xyz):
            raise ValueError(f'{path.name}: nonfinite or invalid coordinates')
    return {'protein_atom_records': len(atom_lines), 'ligand_atom_records': len(ligand_lines)}


def validate_chemistry(path, ligand):
    text = path.read_text()
    if not any(l.strip().lower() == f'data_{ligand}'.lower() for l in text.splitlines()):
        raise ValueError(f'{path.name}: chemical-component identity does not match {ligand}')
    for field in ('_chem_comp_atom.atom_id', '_chem_comp_atom.type_symbol',
                  '_chem_comp_bond.atom_id_1', '_chem_comp_bond.value_order'):
        if field not in text:
            raise ValueError(f'{path.name}: missing chemistry field {field}')


def prepare_data(directory, offline=False):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    manifest = directory / 'manifest.json'
    previous = json.loads(manifest.read_text()) if manifest.exists() else {}
    old_hashes = {s['file']: s['sha256'] for s in previous.get('sources', [])}
    sources, validation = [], []
    for row in DATASETS:
        pdb, ligand = row['pdb'], row['ligand']
        print(f'Dataset: {pdb}, ligand {ligand}', flush=True)
        for name, url in [(f'{pdb}.pdb', row['structure_url']), (f'{ligand}.cif', row['chemistry_url'])]:
            path = directory / name
            download(url, path, offline)
            sha = hashlib.sha256(path.read_bytes()).hexdigest()
            if name in old_hashes and old_hashes[name] != sha:
                raise ValueError(f'{name} changed since the recorded dataset. Use a fresh data directory for a new dataset version.')
            sources.append({'file': name, 'url': url, 'sha256': sha, 'bytes': path.stat().st_size})
        stats = validate_structure(directory / f'{pdb}.pdb', pdb, ligand)
        validate_chemistry(directory / f'{ligand}.cif', ligand)
        validation.append({'pdb': pdb, 'ligand': ligand, 'series': row['series'], **stats})
    content = {'paper': PAPER, 'source': 'RCSB PDB experimental archive',
               'synthetic_fallback': False, 'sources': sources, 'validation': validation}
    temporary = directory / 'manifest.json.tmp'
    temporary.write_text(json.dumps(content, indent=2) + '\n')
    temporary.replace(manifest)
    print('Validated all five crystal structures and five ligand definitions.', flush=True)
    return sources


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--data-dir', type=Path, default=ROOT / 'data')
    args = parser.parse_args()
    try:
        prepare_data(args.data_dir, args.offline)
    except Exception as exc:
        print(f'DATASET ERROR: {exc}', file=sys.stderr)
        sys.exit(1)
