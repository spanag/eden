#!/bin/bash
# Some configuration options are placed here to be accessible from both setup-install-requirements and also regular build-test scripts,
# basically because Homebrew takes too long to install when placed in an uncommon prefix, and also it is not on PATH in the new M1 Mac's.

# LATER make this a supplied parameter if multiarch is needed
BUILD_TARGET_ARCH=$(uname -m)
# LATER check if it's an identified arch. for now assume it matches platform tags on manylinux https://github.com/pypa/manylinux/blob/main/README.rst ; see also https://github.com/pypa/wheel/blob/0.40.0/src/wheel/bdist_wheel.py#L64
if echo $BUILD_TARGET_ARCH | grep -q "^i.86\$" ; then
	BUILD_TARGET_ARCH=i686
fi

# put these here bc manylinux is also decribed like wheel platname
WHEEL_OS_VER=
MANYLINUX_CONTAINER_NEEDED_FOR_WHEEL=

if [ "$BUILD_TARGET_ARCH" = "x86_64" ]; then
	WHEEL_OS_VER="manylinux1"
elif  [ "$BUILD_TARGET_ARCH" = "i686" ]; then
	WHEEL_OS_VER="manylinux1"
elif  [ "$BUILD_TARGET_ARCH" = "arm7l" ]; then
	WHEEL_OS_VER="manylinux2014"
MANYLINUX_CONTAINER_MISSING=1
elif  [ "$BUILD_TARGET_ARCH" = "aarch64" ]; then
	WHEEL_OS_VER="manylinux2014"
	MANYLINUX_CONTAINER_NEEDED_FOR_WHEEL=1
else
	echo "Can't find manylinux version for unknown target arch $BUILD_TARGET_ARCH ! Assume manylinux2014" >&2
	WHEEL_OS_VER="manylinux2014"
fi

WHEEL_TARGET_PLAT="${WHEEL_OS_VER}_${BUILD_TARGET_ARCH}"
WHEEL_PLAT_NAME="linux-${BUILD_TARGET_ARCH}"
WHEEL_PLAT_NAME_FILENAME=$( echo "$WHEEL_PLAT_NAME" | tr .- __ )

MANYLINUX_IMAGE=quay.io/pypa/${WHEEL_OS_VER}_${BUILD_TARGET_ARCH}

# omit homebrew based config until it's needed, refer to mac scripts for that

