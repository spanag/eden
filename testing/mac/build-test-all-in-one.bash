#!/bin/bash
set -e


set PREVPATH="$PATH"

if [ -z "$ARTIFACTS_DIR" ]; then
	ARTIFACTS_DIR=artifacts
	echo "ARTIFACTS_DIR was not defined, set to \"$ARTIFACTS_DIR\""
fi

mkdir -p "$ARTIFACTS_DIR"

# Build for native arch, test, and get wheel
source "$(dirname "${BASH_SOURCE[0]}")/build-wheel.bash"
# REPO_DIR, VERSION, VERSION_PY, WHEEL_PLAT_NAME_FILENAME are set by sourced build script
# echo $REPO_DIR $VERSION_PY $WHEEL_PLAT_NAME $WHEEL_PLAT_NAME_FILENAME
WHEEL_TO_TEST="$REPO_DIR/testing/sandbox/python_package/dist/eden_simulator-$VERSION_PY-py3-none-$WHEEL_PLAT_NAME_FILENAME.whl"

"$(dirname "${BASH_SOURCE[0]}")/run-tests-on-wheel.bash" "$WHEEL_TO_TEST"

cp -f "$WHEEL_TO_TEST" "$ARTIFACTS_DIR"

# Prepare for a new build, reset PATH
set PATH=%PREVPATH%


# done !
