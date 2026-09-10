#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TARGET_BIN="${SCRIPT_DIR}/indelible"

if [ -x "${TARGET_BIN}" ]; then
    # Test if it runs (pipe empty input to avoid interactive cin prompt hang, ignore non-zero exit code)
    CHECK_OUT=$(echo "" | "${TARGET_BIN}" 2>&1 || true)
    if echo "${CHECK_OUT}" | grep -q "INDELible"; then
        echo "INDELible is already installed and functional at: ${TARGET_BIN}"
        exit 0
    fi
fi

echo "Building INDELible v1.03 from evolbioinfo/indelible..."
TMP_BUILD_DIR="$(mktemp -d -t indelible_build_XXXXXX)"
trap 'rm -rf "${TMP_BUILD_DIR}"' EXIT

git clone --depth 1 https://github.com/evolbioinfo/indelible.git "${TMP_BUILD_DIR}/indelible"
cd "${TMP_BUILD_DIR}/indelible/src"

# Compiler detection (g++ or clang++)
CXX="${CXX:-g++}"
if ! command -v "${CXX}" >/dev/null 2>&1; then
    if command -v clang++ >/dev/null 2>&1; then
        CXX="clang++"
    else
        echo "Error: Neither g++ nor clang++ was found in PATH." >&2
        exit 1
    fi
fi

echo "Compiling with ${CXX} -O3..."
${CXX} -O3 -o indelible indelible.cpp -lm

# Install to bin/
cp indelible "${TARGET_BIN}"
chmod +x "${TARGET_BIN}"

echo "INDELible v1.03 successfully installed at: ${TARGET_BIN}"
