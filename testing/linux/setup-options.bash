#!/bin/bash
# Some configuration options are placed here to be accessible from both setup-install-requirements and also regular build-test scripts,
# basically because Homebrew takes too long to install when placed in an uncommon prefix, and also it is not on PATH in the new M1 Mac's.

# LATER make this a supplied parameter if multiarch is needed
BUILD_TARGET_ARCH=$(uname -m)
# LATER check if it's an identified arch. for now assume it mathches platform tags on manylinux https://github.com/pypa/manylinux/blob/main/README.rst

# put these here bc manylinux is also decribed like wheel platname
WHEEL_OS_VER=
if [ "$BUILD_TARGET_ARCH" = "x86_64" ]; then
WHEEL_OS_VER="manylinux1"
elif  [ "$BUILD_TARGET_ARCH" = "i686" ]; then
WHEEL_OS_VER="manylinux1"
elif  [ "$BUILD_TARGET_ARCH" = "arm7l" ]; then
WHEEL_OS_VER="manylinux2014"
elif  [ "$BUILD_TARGET_ARCH" = "aarch64" ]; then
WHEEL_OS_VER="manylinux2014"
else
    echo "Can't find manylinux version for unknown target arch $BUILD_TARGET_ARCH !"
    exit 1 
fi

WHEEL_PLAT_NAME="${WHEEL_OS_VER}-${BUILD_TARGET_ARCH}"
WHEEL_PLAT_NAME_FILENAME=$( echo "$WHEEL_PLAT_NAME" | tr .- __ )

# omit homebrew based config until it's needed, refer to mac scripts for that

