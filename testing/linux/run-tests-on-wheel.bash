#!/bin/bash
set -e

if [ -z "$1" ]; then
	echo "Command needs a command-line argument for wheel to test. Exiting."
	exit 2
fi
if [ -n "$RUN_DIRECT" ]; then

TEST_VENV_PATH=venv-eden-test
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
REPO_DIR="$(dirname "${BASH_SOURCE[0]}")/../../"
# source "$(dirname "${BASH_SOURCE[0]}")/setup-options.bash"

# see https://github.com/pypa/manylinux
# get image name from first platname, right after the last hyphen :D
DOCKER_IMAGE_TEST=$(echo "$1" | rev | cut -d '-' -f 1 | rev | cut -d '.' -f 1 )
# DOCKER_IMAGE_TEST=$WHEEL_OS_VER
echo $DOCKER_IMAGE_TEST

BUILD_PATH=$(realpath $(dirname $1))
# let's put some more junk in build dir...
# ln -s -T /opt/python/cp python3
PIP_INSTALL_EXTRAS=

# if on an ancient platform, use last built wheel versions to save time 
if [ $(echo $DOCKER_IMAGE_TEST | cut -d _ -f 1) == manylinux1 ]; then
	PIP_INSTALL_EXTRAS="numpy==1.21.0"
fi

bash "$(dirname "${BASH_SOURCE[0]}")/docker/sudo_docker.bash" \
run -it --rm --mount type=bind,source=$(realpath ${REPO_DIR}),destination=/repo,readonly --mount "type=bind,source=$BUILD_PATH,destination=/build"  --user $(id -u):$(id -g) \
 -e RUN_DIRECT=1 --workdir /build -e PIP_INSTALL_EXTRAS=$PIP_INSTALL_EXTRAS \
 quay.io/pypa/${DOCKER_IMAGE_TEST} bash -c "set -e; PATH=/build:/opt/python/cp37-cp37m/bin:\$PATH; ln -sfT \$(which cc) gcc; bash /repo/testing/linux/run-tests-on-wheel.bash \"$(basename $1)\""

fi
