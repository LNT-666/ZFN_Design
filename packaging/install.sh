#!/usr/bin/env bash
# ZFN WebSuite — plain-venv install (no docker, no conda).
#
#   ./install.sh            # install into ./.venv and run the smoke test
#   PY=python3.11 ./install.sh
#
# Verified on Ubuntu 22.04 / Python 3.12.3 (lab server ms01).
set -euo pipefail

PY="${PY:-python3}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "[1/4] python: $($PY -V)"

MODE="venv"
if [ ! -d "$ROOT/.venv" ]; then
  # Some slim Debian/Ubuntu images ship python3 without ensurepip.
  if ! "$PY" -m venv "$ROOT/.venv" 2>/dev/null; then
    echo "      python -m venv unavailable (ensurepip missing) -> using --target install"
    MODE="target"
    rm -rf "$ROOT/.venv"
    mkdir -p "$ROOT/.pylibs"
    PYTHONPATH="$ROOT/.pylibs${PYTHONPATH:+:$PYTHONPATH}" \
      "$PY" -m pip install --quiet --target "$ROOT/.pylibs" -e "$ROOT" 2>/dev/null || \
      PYTHONPATH="$ROOT/.pylibs${PYTHONPATH:+:$PYTHONPATH}" \
      "$PY" -m pip install --quiet --target "$ROOT/.pylibs" "$ROOT"
  fi
fi

if [ "$MODE" = "venv" ]; then
  # shellcheck disable=SC1091
  source "$ROOT/.venv/bin/activate"
  PYEXE=python
  echo "[2/4] upgrading pip"
  python -m pip install --quiet --upgrade pip
  echo "[3/4] installing zfn-design (editable) + optional accelerators"
  python -m pip install --quiet -e "$ROOT"
  python -m pip install --quiet "numpy>=1.26" "pyahocorasick>=2.0" || \
    echo "      (accelerators unavailable; the pure-python fallback will be used)"
else
  PYEXE="$PY"
  export PYTHONPATH="$ROOT/.pylibs${PYTHONPATH:+:$PYTHONPATH}"
  echo "[3/4] installed into $ROOT/.pylibs (PYTHONPATH mode)"
fi

echo "[4/4] smoke test"
"$PYEXE" "$ROOT/packaging/smoke_test.py"

echo
echo "done. next:"
if [ "$MODE" = "venv" ]; then
  echo "  source $ROOT/.venv/bin/activate"
else
  echo "  export PYTHONPATH=$ROOT/.pylibs:\$PYTHONPATH"
fi
echo "  zfn-design gui --port 8797 --no-browser --out-root ./out"
