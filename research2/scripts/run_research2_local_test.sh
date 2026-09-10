#!/bin/bash
# Local test runner for research2
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

echo "Running research2 local test runner..."
python3 "${SCRIPT_DIR}/run_research2_local_test.py"
