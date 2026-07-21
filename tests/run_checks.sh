#!/usr/bin/env bash

set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

find "${ROOT}/scripts" "${ROOT}/tests" -name '*.sh' -type f -print0 \
    | sort -z \
    | while IFS= read -r -d '' script; do
        bash -n "${script}"
    done

python3 -m compileall -q "${ROOT}/scripts" "${ROOT}/tests"
python3 -m unittest discover -s "${ROOT}/tests" -p 'test_*.py' -v

echo "All lightweight checks passed."
