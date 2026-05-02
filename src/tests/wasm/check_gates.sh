#!/usr/bin/env bash
# Stage 1 gate check: verify that the WASM-bound translation units compile
# under -DUNIFRAC_WASM=1 -DSKIP_MMAP=1 -DNOGPU=1 with NO HDF5 or lz4 in the
# include path.
#
# This is a syntax-only compile (-fsyntax-only). It mirrors what emcc will
# attempt at Stage 2, but uses plain g++ so we don't need emscripten installed
# to validate gate placement during development.
#
# Run from repo root:   bash src/tests/wasm/check_gates.sh
# Or via make:          make wasm_gate_check
#
# Exit 0 if all WASM-bound TUs parse cleanly without HDF5/lz4. Exit non-zero
# otherwise.

set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "${REPO_ROOT}"

CXX="${CXX:-g++}"
TMPINC="$(mktemp -d)"
trap 'rm -rf "${TMPINC}"' EXIT

# Ensure HDF5 / lz4 cannot leak in via host include path. We construct a
# minimal include set: just the system C++ stdlib + repo src headers. Compiler
# default search paths still include /usr/include, so this script is best-
# effort at simulating the WASM include sandbox; the real validation is the
# emcc compile in Stage 2's wasm_smoke_test.
CFLAGS="-std=c++17 -fsyntax-only -Wall -Wno-unknown-pragmas"
DEFS="-DUNIFRAC_WASM=1 -DSKIP_MMAP=1 -DNOGPU=1"
INC="-I src"

# skbio_alt.cpp pulls <scikit-bio-binaries/...> headers. skbb installs those
# under ${PREFIX}/include/scikit-bio-binaries/ (from src/extern/*.h). If
# SKBB_DIR is set to a sibling skbb checkout, fake the install layout via a
# symlink so `<scikit-bio-binaries/util.h>` resolves. Otherwise skip
# skbio_alt with a notice (Stage 2's emcc build is the authoritative check).
if [[ -n "${SKBB_DIR:-}" && -d "${SKBB_DIR}/src/extern" ]]; then
    ln -s "${SKBB_DIR}/src/extern" "${TMPINC}/scikit-bio-binaries"
    INC="${INC} -I${TMPINC}"
    SKIP_SKBIO_ALT=0
else
    SKIP_SKBIO_ALT=1
fi

# WASM object list per FINDINGS-wasm.md §5. biom.cpp / biom.hpp deliberately
# excluded — the WASM build never reads HDF5 BIOM tables.
WASM_TUS=(
    src/api.cpp
    src/biom_inmem.cpp
    src/biom_subsampled.cpp
    src/skbio_alt.cpp
    src/tree.cpp
    src/unifrac.cpp
    src/unifrac_internal.cpp
    src/unifrac_cmp.cpp
)

fail=0
for tu in "${WASM_TUS[@]}"; do
    if [[ ! -f "${tu}" ]]; then
        echo "SKIP: ${tu} (not present)"
        continue
    fi
    if [[ "${tu}" == "src/skbio_alt.cpp" && "${SKIP_SKBIO_ALT}" == "1" ]]; then
        printf '  SKIP  %-40s (set SKBB_DIR to enable)\n' "${tu}"
        continue
    fi
    printf '  CHECK %-40s ... ' "${tu}"
    # -fsyntax-only means no object output is emitted; -c is intentionally
    # omitted so a stray .o file cannot land in the source tree.
    if ${CXX} ${CFLAGS} ${DEFS} ${INC} "${tu}" 2>"${TMPINC}/err"; then
        echo OK
    else
        echo FAIL
        sed 's/^/      /' "${TMPINC}/err"
        fail=1
    fi
done

if (( fail )); then
    echo
    echo "FAIL: one or more WASM-bound TUs do not compile under -DUNIFRAC_WASM=1"
    exit 1
fi

echo
echo "OK: all WASM-bound TUs parse cleanly under -DUNIFRAC_WASM=1"
