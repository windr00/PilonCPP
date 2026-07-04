#!/bin/bash
# Build a portable piloncpp binary (manylinux_2_28, glibc >= 2.28)
# Compatible with: Ubuntu 18.04+, CentOS 8+, RHEL 8+
#
# Usage: ./build_portable.sh
# Output: piloncpp-manylinux_x86_64

set -euo pipefail

IMAGE="almalinux:8"
HOSTDIR="$(cd "$(dirname "$0")" && pwd)"

echo "=== Building portable piloncpp (manylinux_2_28, glibc >= 2.28) ==="

# Pull image
docker pull "$IMAGE" --quiet 2>/dev/null
echo ""

# Build inside container
docker run --rm \
    -v "${HOSTDIR}":/workspace \
    -w /workspace \
    -e https_proxy="${https_proxy:-}" \
    -e http_proxy="${http_proxy:-}" \
    --platform linux/amd64 \
    "$IMAGE" bash -c '
        set -euo pipefail

        echo "--- Installing build tools ---"
        dnf install -y -q epel-release 2>/dev/null || true
        dnf install -y -q dnf-plugins-core 2>/dev/null || true
        dnf config-manager --set-enabled powertools 2>/dev/null || \
        dnf config-manager --set-enabled crb 2>/dev/null || true
        dnf install -y -q \
            gcc-toolset-13 cmake make zlib-devel \
            2>&1 | tail -3

        source /opt/rh/gcc-toolset-13/enable
        gcc --version | head -1

        echo ""
        echo "--- CMake configure ---"
        rm -rf build_portable && mkdir -p build_portable
        cd build_portable
        cmake .. -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3

        echo ""
        echo "--- Compile ---"
        make -j$(nproc) 2>&1 | tail -5

        echo ""
        echo "--- Strip ---"
        strip piloncpp
        ls -lh piloncpp

        echo ""
        echo "--- Dependencies ---"
        ldd piloncpp 2>&1 || true
'

# Copy binary out of docker volume
if [ -f "${HOSTDIR}/build_portable/piloncpp" ]; then
    cp "${HOSTDIR}/build_portable/piloncpp" "${HOSTDIR}/piloncpp-manylinux_x86_64"
    ls -lh "${HOSTDIR}/piloncpp-manylinux_x86_64"
    echo ""
    echo "=== Done: piloncpp-manylinux_x86_64 ==="
    echo "Compatible with glibc >= 2.28 (Ubuntu 18.04+ / CentOS 8+ / RHEL 8+)"
else
    echo "ERROR: build failed, no binary produced" >&2
    exit 1
fi
