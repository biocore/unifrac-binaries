#!/usr/bin/env bash
# Local mirror of the build-and-test-wasm GitHub Actions job.
#
# Spins up a clean Ubuntu 24.04 container with only the packages the
# CI workflow installs, then runs the same commands the workflow runs.
# Catches the class of bug where a local conda env quietly provides
# headers (lapacke.h, hdf5.h, ...) that the CI runner doesn't have.
#
# Usage:
#   scripts/local_ci_wasm.sh
#
# Prerequisites: docker daemon running, internet access for apt + emsdk
# + GitHub clone of scikit-bio-binaries.
#
# Adapted from scikit-bio-binaries/scripts/local_ci_wasm.sh (removed in
# skbb commit 5d4d183 once their CI was settled). Reintroduced here
# because unifrac's WASM build has a wider dependency surface (skbb
# checkout + LAPACK/BLAS for native expected-value generators) and
# benefits more from a sealed local CI mirror.

set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"

# Pinned skbb ref. Keep in sync with the `ref:` field in
# .github/workflows/main.yml (build-and-test-wasm job).
SKBB_REF="bed8f41002388608c5eeb57d8858dfc76df96ed9"

echo "[local_ci_wasm] mounting ${REPO_ROOT} -> /source (ro) in ubuntu:24.04"
echo "[local_ci_wasm] skbb ref ${SKBB_REF}"

# Mount the repo READ-ONLY at /source and have the container rsync it
# into a private /work directory. This keeps build artifacts container-
# local so the host filesystem is never polluted by docker-as-root files.
docker run --rm -t \
  --network=host \
  -v "${REPO_ROOT}:/source:ro" \
  -e DEBIAN_FRONTEND=noninteractive \
  -e SKBB_REF="${SKBB_REF}" \
  ubuntu:24.04 \
  bash -ec '
    set -eu
    export EMSDK_QUIET=1

    echo "=== container identity ==="
    id
    uname -a
    cat /etc/os-release | head -3

    # Copy the repo into a private writable location, excluding
    # build artifacts (force a fresh build inside the container so
    # we test the genuine bootstrap, not the host-cached state).
    mkdir -p /work/unifrac-binaries
    cp -a /source/. /work/unifrac-binaries/
    rm -rf /work/unifrac-binaries/src/.wasm-skbb-include
    rm -f  /work/unifrac-binaries/src/libssu_wasm.a
    rm -f  /work/unifrac-binaries/src/*.wasm.o
    rm -f  /work/unifrac-binaries/src/*_wasm.js /work/unifrac-binaries/src/*_wasm.wasm
    rm -f  /work/unifrac-binaries/src/unifrac_accapi_cpu.cpp
    rm -f  /work/unifrac-binaries/src/unifrac_task_noclass_cpu.cpp

    # Step 5 in main.yml — apt installs. Keep this list in sync with
    # the apt-get install line in .github/workflows/main.yml
    # (build-and-test-wasm job); if they drift, this harness no
    # longer mirrors CI.
    #
    # Required packages:
    #   ca-certificates, curl, git, python3, xz-utils — checkout/emsdk
    #   g++, make                                     — native generators
    #   libopenblas-dev, libblas-dev, liblapacke-dev — native skbb +
    #                                                  unifrac native
    #                                                  generators (Stages
    #                                                  5/6/7).
    apt-get update -qq
    apt-get install -y --no-install-recommends \
        ca-certificates curl git python3 xz-utils \
        g++ make \
        libopenblas-dev libblas-dev liblapacke-dev

    # Steps 2/3 in main.yml — checkout sibling skbb at the pinned
    # commit + setup-emsdk + setup-node.
    cd /work
    git clone --depth 1 https://github.com/scikit-bio/scikit-bio-binaries.git scikit-bio-binaries
    ( cd scikit-bio-binaries && git fetch --depth=1 origin "${SKBB_REF}" && git checkout "${SKBB_REF}" )

    if [ ! -d /opt/emsdk ]; then
      git clone --depth 1 https://github.com/emscripten-core/emsdk.git /opt/emsdk
    fi
    /opt/emsdk/emsdk install 3.1.71
    /opt/emsdk/emsdk activate 3.1.71
    # shellcheck source=/dev/null
    . /opt/emsdk/emsdk_env.sh

    # Workaround for an emsdk + GNU make interaction:
    # emsdk_env.sh prepends /opt/emsdk to PATH, and that directory
    # contains a `node` SUBDIRECTORY (the toolchain root for the
    # bundled node). When make does a PATH lookup for "node" via its
    # direct-execve fast path, it finds that directory FIRST and tries
    # to exec it, which returns EACCES. bash filters non-regular files
    # during PATH lookup and avoids this; CI uses actions/setup-node
    # which prepends a real node-bin directory. We mirror that.
    export PATH="$(dirname "${EMSDK_NODE}"):${PATH}"
    node --version

    # Steps 8/9/10 — build skbb WASM, then unifrac WASM, then run tests.
    ( cd scikit-bio-binaries && make wasm )
    ( cd unifrac-binaries     && make wasm SKBB_DIR=/work/scikit-bio-binaries )
    ( cd unifrac-binaries     && make wasm_test SKBB_DIR=/work/scikit-bio-binaries )

    # Sanity-check the artifact
    ls -la /work/unifrac-binaries/src/libssu_wasm.a
    emar t   /work/unifrac-binaries/src/libssu_wasm.a | head -20

    echo
    echo "[local_ci_wasm] PASSED — local CI mirror is green"
  '
