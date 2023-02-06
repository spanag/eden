#!/bin/bash
set -e

if [ -z "$1" ]; then
	echo "Command needs a command-line argument for wheel to test. Exiting."
	exit 2
fi

TEST_VENV_PATH=venv-eden-test

# Use a clean venv
rm -rf "$TEST_VENV_PATH"
python3 -m venv "$TEST_VENV_PATH"
source  "$TEST_VENV_PATH/bin/activate"

python3 -m pip uninstall -y eden-simulator 
python3 -m pip install wheel
python3 -m pip install "$1"

source "$(dirname "${BASH_SOURCE[0]}")/run-basic-tests.bash"

