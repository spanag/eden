#!/bin/bash
set -e

source "$(dirname "${BASH_SOURCE[0]}")/default-repo-path.bash"

# more elegant LATER
source "$(dirname "${BASH_SOURCE[0]}")/setup-options.bash"

source "$(dirname "${BASH_SOURCE[0]}")/get-version.bash"
VERSION_PY=$VERSION

make -f "${REPO_DIR}/testing/linux/docker/Makefile" docker_build_env

# rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"

# TODO mount source tree readonly
bash "$(dirname "${BASH_SOURCE[0]}")/docker/sudo_docker.bash" \
run -it --rm --mount type=bind,source=${REPO_DIR},destination=/repo --mount "type=bind,source=$(realpath $BUILD_DIR),destination=/build"  --user $(id -u):$(id -g) \
-e OUT_DIR=/build -e TARGETS=wheel -e BUILD=release -e BUILD_STAMP="$VERSION" -e WHEEL_TARGET_PLAT="$WHEEL_PLAT_NAME_FILENAME" -e WHEEL_PLAT_NAME_FILENAME="$WHEEL_PLAT_NAME_FILENAME" -e EXTRA_WHEEL_PACKAGE_TAGS="--plat-name $WHEEL_PLAT_NAME" -e CFLAGS_extra="-static" --workdir /repo \
 eden-build:latest bash /repo/testing/linux/docker/build_on_docker.bash

# make -j$(sysctl -n hw.logicalcpu) eden wheel BUILD=release BUILD_STAMP="$VERSION" WHEEL_VERSION="$VERSION_PY" WHEEL_PLAT="$WHEEL_PLAT_NAME" WHEEL_TARGET_PLAT="$WHEEL_PLAT_NAME" WHEEL_PLAT_NAME_FILENAME="$WHEEL_PLAT_NAME_FILENAME" EXTRA_WHEEL_PACKAGE_TAGS="--plat-name $WHEEL_PLAT_NAME"

