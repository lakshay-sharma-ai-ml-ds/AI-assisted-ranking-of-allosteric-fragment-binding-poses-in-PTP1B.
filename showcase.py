"""Refresh README results from a completed run, without rerunning the science."""
import csv
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
START, END = '<!-- RESULTS:START -->', '<!-- RESULTS:END -->'
PDBS = {'8G65', '8G68', '8G69', '8G67', '8G6A'}
METHODS = ('Vina', 'ML', 'Shuffled-label ML')


def result_section(directory):
    directory = Path(directory)
    required = ('metrics.csv', 'summary.json', 'manifest.json', 'REPORT.md', 'ranking.png', 'poses.png')
    missing = [name for name in required if not (directory / name).is_file()]
    if missing:
        raise ValueError('No completed result set: missing ' + ', '.join(missing))
    manifest = json.loads((directory / 'manifest.json').read_text())
    if 'known-pocket holo redocking' not in manifest.get('protocol', ''):
        raise ValueError('Results do not identify the expected experiment protocol')
    source_files = {s['file'] for s in manifest.get('sources', [])}
    if not {p + '.pdb' for p in PDBS}.issubset(source_files):
        raise ValueError('Run manifest is missing experimental structure provenance')
    with (directory / 'metrics.csv').open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    keyed = {(r['fragment'], r['method']): r for r in rows}
    if len(rows) != 15 or len(keyed) != 15 or set(keyed) != {(p, m) for p in PDBS for m in METHODS}:
        raise ValueError('Expected all five fragments and all three comparison methods')
    for r in rows:
        for k in ('top1_rmsd_A', 'oracle_rmsd_A'):
            v = float(r[k])
            if not math.isfinite(v) or v < 0:
                raise ValueError('Invalid RMSD in metrics')
        if int(r['n_poses']) < 1 or r['top1_success'] not in ('0', '1') or r['top3_success'] not in ('0', '1'):
            raise ValueError('Invalid pose count or success label')
    summary = json.loads((directory / 'summary.json').read_text())
    delta = summary['paired_ML_minus_Vina_top1_rmsd_A']
    difference = sum(float(keyed[p, 'ML']['top1_rmsd_A']) - float(keyed[p, 'Vina']['top1_rmsd_A'])
                     for p in PDBS) / len(PDBS)
    if not math.isclose(difference, float(delta['mean_difference']), abs_tol=1e-8):
        raise ValueError('Metrics and summary disagree; use files from the same run')
    interval = list(map(float, delta['percentile_95_CI']))
    if len(interval) != 2 or not all(math.isfinite(v) for v in interval) or interval[0] > interval[1]:
        raise ValueError('Invalid comparison interval')
    lines = ['### Measured pose recovery', '',
             '| Fragment | Candidates | Vina top-1 RMSD (Å) | ML top-1 RMSD (Å) | Best candidate RMSD (Å) |',
             '|---|---:|---:|---:|---:|']
    for pdb in sorted(PDBS):
        a, b = keyed[pdb, 'Vina'], keyed[pdb, 'ML']
        lines.append(f"| {pdb} | {a['n_poses']} | {float(a['top1_rmsd_A']):.3f} | {float(b['top1_rmsd_A']):.3f} | {float(a['oracle_rmsd_A']):.3f} |")
    lines += ['', '| Method | Top-1 recovery ≤ 2 Å | Top-3 recovery ≤ 2 Å |', '|---|---:|---:|']
    for method in METHODS:
        selected = [keyed[p, method] for p in PDBS]
        top1 = sum(int(r['top1_success']) for r in selected)
        top3 = sum(int(r['top3_success']) for r in selected)
        lines.append(f'| {method} | {top1}/5 | {top3}/5 |')
    lines += ['', f'**Mean paired ML − Vina top-1 RMSD: {difference:+.3f} Å.** Negative favors ML.',
              f'Descriptive fragment-bootstrap 95% interval: [{interval[0]:.3f}, {interval[1]:.3f}] Å.',
              'Five related fragments across two chemical series; these results are exploratory.', '',
              '![Held-out pose-ranking comparison](results/ranking.png)', '',
              '*Vina versus held-out-series ML, with candidate scores and crystal-pose error.*', '',
              '![Crystal and ML-ranked fragment poses](results/poses.png)', '',
              '*Experimental crystal poses and ML top-ranked candidates in the same receptor coordinate frame.*', '',
              '[Full analysis](results/REPORT.md) · [Metrics](results/metrics.csv) · [Run provenance](results/manifest.json)']
    return '\n'.join(lines)


def update_readme(readme=ROOT / 'README.md', results=ROOT / 'results'):
    readme = Path(readme)
    original = readme.read_text()
    if original.count(START) != 1 or original.count(END) != 1 or original.index(START) >= original.index(END):
        raise ValueError('README must contain one ordered pair of results markers')
    section = result_section(results)
    before, rest = original.split(START)
    _, after = rest.split(END)
    temporary = readme.with_name(readme.name + '.tmp')
    temporary.write_text(before + START + '\n' + section + '\n' + END + after)
    temporary.replace(readme)


if __name__ == '__main__':
    try:
        update_readme()
        print('README updated with measured results and plot links.')
    except (ValueError, KeyError, OSError) as exc:
        print(f'README RESULTS ERROR: {exc}', file=sys.stderr)
        sys.exit(1)
