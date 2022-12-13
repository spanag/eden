#!/bin/bash
set -e

source "$(dirname "${BASH_SOURCE[0]}")/default-repo-path.bash"

# more elegant LATER
source "$(dirname "${BASH_SOURCE[0]}")/setpath.bash"

source "$(dirname "${BASH_SOURCE[0]}")/get-version.bash"
VERSION_PY=$VERSION

cd $REPO_DIR

TARGET_ARCH="$(uname -m)"

WHEEL_MACOS_VER=
if [ "$TARGET_ARCH" = "x86_64" ]; then
WHEEL_MACOS_VER="-10.6"
elif  [ "$TARGET_ARCH" = "arm64" ]; then
WHEEL_MACOS_VER="-11.0"
fi

WHEEL_PLAT_NAME="macosx${WHEEL_MACOS_VER}-$(uname -m)"
WHEEL_PLAT_NAME_FILENAME=$( echo "$WHEEL_PLAT_NAME" | tr .- __ )

make clean
make -j$(sysctl -n hw.logicalcpu) eden wheel BUILD=release BUILD_STAMP="$VERSION" WHEEL_VERSION="$VERSION_PY" WHEEL_PLAT="$WHEEL_PLAT_NAME" WHEEL_TARGET_PLAT="$WHEEL_PLAT_NAME" WHEEL_PLAT_NAME_FILENAME="$WHEEL_PLAT_NAME_FILENAME" EXTRA_WHEEL_PACKAGE_TAGS="--plat-name $WHEEL_PLAT_NAME"

cd -

