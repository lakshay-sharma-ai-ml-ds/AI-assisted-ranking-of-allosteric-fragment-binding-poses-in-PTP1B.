"""Launcher contracts and README generation, using explicitly artificial fixtures."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from showcase import PDBS, METHODS, update_readme


def launcher(tmp_path, arguments=(), failure=False):
    python = tmp_path / 'interpreter'
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$*" >> "$TRACE"\n'
                      'if [[ "$1" == "pipeline.py" && "${FAIL_PIPELINE:-}" == "yes" ]]; then exit 9; fi\n')
    python.chmod(0o755)
    trace = tmp_path / 'trace'
    env = dict(os.environ, PYTHON=str(python), TRACE=str(trace), FAIL_PIPELINE='yes' if failure else 'no')
    result = subprocess.run(['bash', str(ROOT / 'run.sh'), *arguments], env=env, capture_output=True, text=True)
    return result, trace.read_text().splitlines()


def test_launcher_uses_supplied_interpreter_and_forwards_args(tmp_path):
    result, calls = launcher(tmp_path, ['--offline', '--cpu', '4'])
    assert result.returncode == 0
    assert calls == ['-', 'pipeline.py --offline --cpu 4', 'showcase.py']


def test_launcher_failed_pipeline_does_not_publish_results(tmp_path):
    result, calls = launcher(tmp_path, failure=True)
    assert result.returncode == 9
    assert 'showcase.py' not in calls


def test_download_only_does_not_check_scientific_dependencies(tmp_path):
    result, calls = launcher(tmp_path, ['--download-only', '--offline'])
    assert result.returncode == 0
    assert calls == ['dataset.py --offline']


def test_launcher_does_not_install_or_create_environment():
    script = (ROOT / 'run.sh').read_text()
    assert 'pip install' not in script
    assert 'conda' not in script.lower()
    assert 'miniforge' not in script.lower()
    assert '-m venv' not in script


def result_fixture(directory):
    # This fixture exists only inside the temporary test directory.
    directory.mkdir()
    rows = []
    for pdb in sorted(PDBS):
        for method in METHODS:
            rmsd = 1.0 if method == 'ML' else 2.0
            rows.append({'fragment': pdb, 'method': method, 'n_poses': 4,
                         'top1_rmsd_A': rmsd, 'oracle_rmsd_A': 1,
                         'top1_success': 1, 'top3_success': 1})
    with (directory / 'metrics.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (directory / 'manifest.json').write_text(json.dumps({'protocol': 'known-pocket holo redocking',
                            'sources': [{'file': p + '.pdb'} for p in PDBS]}))
    (directory / 'summary.json').write_text(json.dumps({'paired_ML_minus_Vina_top1_rmsd_A':
                                {'mean_difference': -1, 'percentile_95_CI': [-1, -1]}}))
    for filename in ('REPORT.md', 'ranking.png', 'poses.png'):
        (directory / filename).write_text('TEST FIXTURE ONLY')


def test_results_are_embedded_without_changing_other_readme_sections(tmp_path):
    results = tmp_path / 'results'
    result_fixture(results)
    readme = tmp_path / 'README.md'
    readme.write_text('Intro\n<!-- RESULTS:START -->\nPending\n<!-- RESULTS:END -->\nArchitecture\n')
    update_readme(readme, results)
    assert readme.read_text().startswith('Intro\n')
    assert readme.read_text().endswith('\nArchitecture\n')
    assert '-1.000 Å' in readme.read_text()
    assert '](results/ranking.png)' in readme.read_text()
    once = readme.read_text()
    update_readme(readme, results)
    assert readme.read_text() == once


def test_incomplete_results_leave_readme_untouched(tmp_path):
    results = tmp_path / 'results'
    results.mkdir()
    readme = tmp_path / 'README.md'
    original = '<!-- RESULTS:START -->Pending<!-- RESULTS:END -->'
    readme.write_text(original)
    with pytest.raises(ValueError, match='No completed'):
        update_readme(readme, results)
    assert readme.read_text() == original


def test_mixed_run_metrics_and_summary_rejected(tmp_path):
    results = tmp_path / 'results'
    result_fixture(results)
    (results / 'summary.json').write_text(json.dumps({'paired_ML_minus_Vina_top1_rmsd_A':
                                {'mean_difference': 4, 'percentile_95_CI': [3, 5]}}))
    readme = tmp_path / 'README.md'
    readme.write_text('<!-- RESULTS:START -->Pending<!-- RESULTS:END -->')
    with pytest.raises(ValueError, match='disagree'):
        update_readme(readme, results)
