#!/usr/bin/env bash

set -euo pipefail

readonly SDSL_VERSION="2.1.1"
readonly SDSL_COMMIT="0546faf0552142f06ff4b201b671a5769dd007ad"
readonly SDSL_SOURCE_DIR="${PIXI_PROJECT_ROOT}/.pixi/source/sdsl-lite-${SDSL_VERSION}"
readonly SDSL_BUILD_DIR="${PIXI_PROJECT_ROOT}/.pixi/build/sdsl-lite-${SDSL_VERSION}"
readonly SDSL_LIBRARY="${CONDA_PREFIX}/lib/libsdsl.a"

if [[ -f "${SDSL_LIBRARY}" ]]; then
    exit 0
fi

mkdir -p "$(dirname "${SDSL_SOURCE_DIR}")" "$(dirname "${SDSL_BUILD_DIR}")"

if [[ ! -d "${SDSL_SOURCE_DIR}/.git" ]]; then
    git clone \
        --branch "v${SDSL_VERSION}" \
        --depth 1 \
        --recurse-submodules \
        --shallow-submodules \
        https://github.com/simongog/sdsl-lite.git \
        "${SDSL_SOURCE_DIR}"
fi

actual_commit="$(git -C "${SDSL_SOURCE_DIR}" rev-parse HEAD)"
if [[ "${actual_commit}" != "${SDSL_COMMIT}" ]]; then
    echo "Expected SDSL commit ${SDSL_COMMIT}, found ${actual_commit}." >&2
    echo "Remove ${SDSL_SOURCE_DIR} and rerun this task." >&2
    exit 1
fi

git -C "${SDSL_SOURCE_DIR}" submodule update --init --recursive

readonly SDSL_PATCH="${PIXI_PROJECT_ROOT}/scripts/sdsl-2.1.1-apple-clang.patch"
if git -C "${SDSL_SOURCE_DIR}" apply --check "${SDSL_PATCH}" 2>/dev/null; then
    git -C "${SDSL_SOURCE_DIR}" apply "${SDSL_PATCH}"
elif ! git -C "${SDSL_SOURCE_DIR}" apply --reverse --check "${SDSL_PATCH}" 2>/dev/null; then
    echo "SDSL compatibility patch cannot be applied cleanly." >&2
    exit 1
fi

BUILD_PORTABLE=1 cmake \
    -S "${SDSL_SOURCE_DIR}" \
    -B "${SDSL_BUILD_DIR}" \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${CONDA_PREFIX}"

cmake --build "${SDSL_BUILD_DIR}" \
    --target sdsl divsufsort divsufsort64 \
    --parallel
cmake --install "${SDSL_BUILD_DIR}"
