#!/bin/bash
set -e


set PREVPATH="$PATH"

if [ -z "$ARTIFACTS_DIR" ]; then
	ARTIFACTS_DIR=artifacts
	echo "ARTIFACTS_DIR was not defined, set to \"$ARTIFACTS_DIR\""
fi
if [ -z "$BUILD_DIR" ]; then
	BUILD_DIR=build
	echo "BUILD_DIR was not defined, set to \"$BUILD_DIR\""
fi

mkdir -p "$ARTIFACTS_DIR"

# Build for native arch, test, and get wheel
source "$(dirname "${BASH_SOURCE[0]}")/build-wheel.bash"
# REPO_DIR, VERSION, VERSION_PY, WHEEL_PLAT_NAME_FILENAME are set by sourced build script
# echo $REPO_DIR $VERSION_PY $WHEEL_PLAT_NAME $WHEEL_PLAT_NAME_FILENAME
WHEEL_TO_TEST=$(find "$BUILD_DIR/bin" -type f -name "eden_simulator-$VERSION_PY-py3-none-*.whl")
echo $WHEEL_TO_TEST

"$(dirname "${BASH_SOURCE[0]}")/run-tests-on-wheel.bash" "$WHEEL_TO_TEST"
echo "Tests passed"
cp -f "$WHEEL_TO_TEST" "$ARTIFACTS_DIR"
echo "Artifact ready on $ARTIFACTS_DIR"

# done !
