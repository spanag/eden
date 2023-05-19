#!/bin/bash
set -e

if [ -z "$1" ]; then
	echo "Command needs a command-line argument for wheel to test. Exiting."
	exit 2
fi

# https://stackoverflow.com/questions/5265817/how-to-get-full-path-of-a-file
BUILD_PATH=$(python3 -c "import os; print(os.path.abspath(\"$(dirname $1)\"))")

if [ -n "$RUN_DIRECT" ]; then

TEST_VENV_PATH=$BUILD_PATH/venv-eden-test
# gcc
# python3
# Use a clean venv
rm -rf "$TEST_VENV_PATH"
# TEST_VENV_PATH=$(mktemp -d)
python3 -m venv "$TEST_VENV_PATH"
source  "$TEST_VENV_PATH/bin/activate"
python3 -m pip install -U pip
python3 -m pip uninstall -y eden-simulator 

python3 -m pip install $PIP_INSTALL_EXTRAS wheel
python3 -m pip install "$1"

source "$(dirname "${BASH_SOURCE[0]}")/run-basic-tests.bash"

else

# run on Docker
source "$(dirname "${BASH_SOURCE[0]}")/setup-options.bash"
source "$(dirname "${BASH_SOURCE[0]}")/docker/setup-options.bash"

if [ -n "$MANYLINUX_CONTAINER_MISSING" ]; then
	# just run it directly after all :c
	RUN_DIRECT=1 ${BASH_SOURCE[0]} "$@"
	exit $? # just in case -e is not set
fi

REPO_DIR="$(dirname "${BASH_SOURCE[0]}")/../../"

PIP_INSTALL_EXTRAS=

# if on an ancient platform, use last built wheel versions to save time 
if [ $WHEEL_OS_VER == manylinux1 ] || [ $WHEEL_OS_VER == manylinux2010 ]; then
	PIP_INSTALL_EXTRAS="numpy==1.21.0"
fi

bash "$(dirname "${BASH_SOURCE[0]}")/docker/sudo_docker.bash" \
run -it --rm --mount type=bind,source=$(realpath ${REPO_DIR}),destination=/repo,readonly --mount "type=bind,source=$BUILD_PATH,destination=/build"  --user $DOCKER_USER_OR_ROOT \
 -e RUN_DIRECT=1 --workdir /build -e PIP_INSTALL_EXTRAS=$PIP_INSTALL_EXTRAS \
 $MANYLINUX_IMAGE bash -c "set -e; PATH=/build:/opt/python/cp37-cp37m/bin:\$PATH; ln -sfT \$(which cc) gcc; bash /repo/testing/linux/run-tests-on-wheel.bash \"$(basename $1)\""

fi
