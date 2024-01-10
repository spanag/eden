#!/bin/bash
set -e

source "$(dirname "${BASH_SOURCE[0]}")/default-repo-path.bash"

# more elegant LATER
source "$(dirname "${BASH_SOURCE[0]}")/setup-options.bash"
source "$(dirname "${BASH_SOURCE[0]}")/docker/setup-options.bash"

source "$(dirname "${BASH_SOURCE[0]}")/get-version.bash"

WHEEL_VERSION=$VERSION

rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"

#TODO add an option to build *without* using Docker
# and merge with the 'build docs' use case

BUILD_ENVVAR_DOCKER_EXTRA=
if [ -n "$MANYLINUX_CONTAINER_NEEDED_FOR_WHEEL" ]; then
	BUILD_ENV_IMAGE=eden-manylinux:latest
	BUILD_ENV_MAKEFILE_TARGET=docker_build_manylinux

	export MANYLINUX_IMAGE
	export MANYLINUX_VERSION=$WHEEL_OS_VER
else
	BUILD_ENV_IMAGE=eden-build:latest
	BUILD_ENV_MAKEFILE_TARGET=docker_build_env

	BUILD_ENVVAR_DOCKER_EXTRA="$BUILD_ENVVAR_DOCKER_EXTRA -e CFLAGS_extra=\"-static\""
fi

make -f "${REPO_DIR}/testing/linux/docker/Makefile" $BUILD_ENV_MAKEFILE_TARGET

# TODO mount source tree readonly - needs to remove use of sandbox...
$SUDO_DOCKER run -it --rm \
--mount type=bind,source=${REPO_DIR},destination=/repo --mount "type=bind,source=$(realpath $BUILD_DIR),destination=/build" \
--user $DOCKER_USER_OR_ROOT \
-e OUT_DIR=/build -e TARGETS=wheel -e BUILD=release -e BUILD_STAMP="$VERSION" -e WHEEL_VERSION="$WHEEL_VERSION" -e WHEEL_TARGET_PLAT="$WHEEL_TARGET_PLAT" -e WHEEL_PLAT_NAME_FILENAME="$WHEEL_PLAT_NAME_FILENAME" $BUILD_ENVVAR_DOCKER_EXTRA -e EXTRA_WHEEL_PACKAGE_TAGS="--plat-name $WHEEL_PLAT_NAME" --workdir /repo \
$BUILD_ENV_IMAGE bash /repo/testing/linux/docker/build_on_docker.bash

# make -j$(sysctl -n hw.logicalcpu) eden wheel BUILD=release BUILD_STAMP="$VERSION" WHEEL_VERSION="$VERSION" WHEEL_PLAT="$WHEEL_PLAT_NAME" WHEEL_TARGET_PLAT="$WHEEL_PLAT_NAME" WHEEL_PLAT_NAME_FILENAME="$WHEEL_PLAT_NAME_FILENAME" EXTRA_WHEEL_PACKAGE_TAGS="--plat-name $WHEEL_PLAT_NAME"

