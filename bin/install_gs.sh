#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TARGET_BIN="${SCRIPT_DIR}/gs2"

if [ -x "${TARGET_BIN}" ] && "${TARGET_BIN}" -v >/dev/null 2>&1; then
    echo "gs2 is already installed and functional at: ${TARGET_BIN}"
    exit 0
fi

echo "Building gs2 from MotomuMatsui/gs..."
TMP_BUILD_DIR="$(mktemp -d -t gs_build_XXXXXX)"
trap 'rm -rf "${TMP_BUILD_DIR}"' EXIT

git clone --depth 1 https://github.com/MotomuMatsui/gs "${TMP_BUILD_DIR}/gs"
cd "${TMP_BUILD_DIR}/gs"

mkdir -p lib
tar xzf lapack-3.12.1.tar.gz
cp lapack-3.12.1/LAPACKE/include/*.h lib/
cp lapack-3.12.1/CBLAS/include/*.h lib/ 2>/dev/null || true

cat << 'EOF' > lib/lapacke_mangling.h
#ifndef LAPACKE_MANGLING_H
#define LAPACKE_MANGLING_H
#ifndef LAPACK_GLOBAL
#define LAPACK_GLOBAL(lcname,UCNAME)  lcname##_
#define LAPACK_GLOBAL_(lcname,UCNAME) lcname##_
#endif
#endif
EOF

cat << 'EOF' > lib/cblas_mangling.h
#ifndef CBLAS_MANGLING_H
#define CBLAS_MANGLING_H
#ifndef CBLAS_GLOBAL
#define CBLAS_GLOBAL(lcname,UCNAME)  lcname##_
#define CBLAS_GLOBAL_(lcname,UCNAME) lcname##_
#endif
#endif
EOF

# Detect OS
OS="$(uname -s)"
ENV_LIB="${CONDA_PREFIX:-${HOME}/.micromamba/envs/phylomethod_env}/lib"

if [ "${OS}" = "Darwin" ]; then
    clang++ -O3 -std=c++14 -march=native -fno-exceptions -funroll-loops -Wall \
        -Ilib -Isrc \
        src/messages.cpp src/eigen.cpp src/transitivity.cpp src/format.cpp src/mmseqs.cpp src/blastn.cpp \
        src/gs_functions.cpp src/sc_functions.cpp src/sc.cpp src/ep.cpp src/gs.cpp src/main.cpp \
        -L"${ENV_LIB}" -Wl,-rpath,"${ENV_LIB}" \
        -llapacke -lopenblas -lm \
        -o "${TARGET_BIN}"
else
    # Linux (LSF supercomputer / HPC)
    g++ -O3 -std=c++14 -march=native -fno-exceptions -funroll-loops -fopenmp -Wall \
        -Ilib -Isrc \
        src/messages.cpp src/eigen.cpp src/transitivity.cpp src/format.cpp src/mmseqs.cpp src/blastn.cpp \
        src/gs_functions.cpp src/sc_functions.cpp src/sc.cpp src/ep.cpp src/gs.cpp src/main.cpp \
        -L"${ENV_LIB}" -Wl,-rpath,"${ENV_LIB}" \
        -llapacke -lopenblas -lm -lgomp \
        -o "${TARGET_BIN}" || \
    g++ -O3 -std=c++14 -fno-exceptions -Wall \
        -Ilib -Isrc \
        src/messages.cpp src/eigen.cpp src/transitivity.cpp src/format.cpp src/mmseqs.cpp src/blastn.cpp \
        src/gs_functions.cpp src/sc_functions.cpp src/sc.cpp src/ep.cpp src/gs.cpp src/main.cpp \
        -L"${ENV_LIB}" -Wl,-rpath,"${ENV_LIB}" \
        -llapacke -lblas -llapack -lm \
        -o "${TARGET_BIN}"
fi

chmod +x "${TARGET_BIN}"
echo "gs2 successfully built and installed at: ${TARGET_BIN}"
"${TARGET_BIN}" -v
