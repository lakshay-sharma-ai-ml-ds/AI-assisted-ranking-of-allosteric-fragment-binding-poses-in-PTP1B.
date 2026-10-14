#!/usr/bin/env bash
# Use the caller's Python environment. Never create one or install packages.
set -euo pipefail
cd -- "$(dirname -- "$0")"

if [[ -z "${PYTHON:-}" ]]; then
    if command -v python >/dev/null 2>&1; then
        PYTHON=python
    elif command -v python3 >/dev/null 2>&1; then
        PYTHON=python3
    else
        echo 'Python is missing. Activate your prepared Python environment first.' >&2
        exit 2
    fi
fi

if [[ "${1:-}" == '--download-only' ]]; then
    shift
    exec "$PYTHON" dataset.py "$@"
fi

export MPLCONFIGDIR="$PWD/.cache/matplotlib"
mkdir -p "$MPLCONFIGDIR"

"$PYTHON" - <<'PY'
import importlib
import shutil
import sys
missing = []
for module in ('numpy', 'scipy', 'sklearn', 'matplotlib', 'rdkit', 'Bio', 'meeko', 'gemmi'):
    try:
        importlib.import_module(module)
    except Exception as exc:
        missing.append(f'{module}: {exc}')
for binary in ('vina', 'obabel'):
    if not shutil.which(binary):
        missing.append(f'{binary}: executable not on PATH')
if missing:
    print('Dependency check failed:\n  ' + '\n  '.join(missing), file=sys.stderr)
    print('Install Python packages from requirements.txt in your environment.', file=sys.stderr)
    print('Install Vina and Open Babel separately and add their executables to PATH.', file=sys.stderr)
    sys.exit(2)
print('Using Python:', sys.executable)
PY

if [[ "${1:-}" == '--check-only' ]]; then exit 0; fi
if [[ "${1:-}" == '--test' ]]; then
    shift
    exec "$PYTHON" -m pytest -q tests "$@"
fi

"$PYTHON" pipeline.py "$@"
"$PYTHON" showcase.py
