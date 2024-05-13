#!/bin/bash
set -e

set PREVPATH="$PATH"
# no, $"\n" won't work: https://unix.stackexchange.com/questions/214858/new-line-in-bash-variables
LF="
"
# set -v
# RUN_DIRECT=1

DEBUG_SPHINX= # for debugging
if [ -n "$BUILD_IN_PLACE" ]; then
	echo "BUILD_IN_PLACE set"
	if [ -z "$ARTIFACTS_DIR" ]; then
		ARTIFACTS_DIR=.
		echo "ARTIFACTS_DIR was not defined, set to \"$ARTIFACTS_DIR\""
	fi
	if [ -z "$BUILD_DIR" ]; then
		BUILD_DIR=.
		echo "BUILD_DIR was not defined, set to \"$BUILD_DIR\""
	fi
else
	if [ -z "$ARTIFACTS_DIR" ]; then
		ARTIFACTS_DIR=artifacts_docs
		echo "ARTIFACTS_DIR was not defined, set to \"$ARTIFACTS_DIR\""
	fi
	if [ -z "$BUILD_DIR" ]; then
		BUILD_DIR=build_docs
		echo "BUILD_DIR was not defined, set to \"$BUILD_DIR\""
	fi
fi

# https://stackoverflow.com/questions/5265817/how-to-get-full-path-of-a-file
BUILD_DIR=$(python3 -c "import os; print(os.path.abspath(\"$BUILD_DIR\"))")

mkdir -p "$ARTIFACTS_DIR" "$BUILD_DIR"


# now set up the config...
source "$(dirname "${BASH_SOURCE[0]}")/default-repo-path.bash"

# more elegant LATER
source "$(dirname "${BASH_SOURCE[0]}")/setup-options.bash"
source "$(dirname "${BASH_SOURCE[0]}")/get-version.bash"
# Build up the line delimited build env options (LATER use \0 if really needed), to use either with env or with docker:
# echo $BUILD_DIR
# exit 0
ENV+=OUT_DIR="$BUILD_DIR"$LF
ENV+=BUILD_STAMP="$VERSION"$LF
ENV+=BUILD="release"$LF

ENV+=TARGETS="eden wheel"$LF
ENV+=WHEEL_VERSION="$VERSION"$LF
ENV+=WHEEL_DONT_REPAIR="true"$LF # instead of qualifying a plat_name

if [ -n "$RUN_DIRECT" ]; then
	
	# skip the venv if deep in thought ...
	if [ -z "$DEEP_IN_THOUGHT" ]; then

		# Make the venv for building the wheel AND using it to run the docs notebooks
		PIP_INSTALL_BUILD_EDEN="wheel auditwheel"
		PIP_INSTALL_BUILD_DOCS_EXTRA= # the rest are in requirements.txt

		TEST_VENV_PATH=$BUILD_DIR/venv-eden-test #TODO rename?

		# flag for debugging sphinx...
		if [ -z "$DEBUG_SPHINX" ]; then
			# Use a clean venv
			rm -rf "$TEST_VENV_PATH"
			# TEST_VENV_PATH=$(mktemp -d)
			python3 -m venv "$TEST_VENV_PATH"
		fi

		source "$TEST_VENV_PATH/bin/activate"
		pip list
		# bash
		if [ -z "$DEBUG_SPHINX" ]; then
			python3 -m pip install -U pip 
			python3 -m pip install $PIP_INSTALL_BUILD_EDEN

			# Build for native arch, get wheel, and run docs
			# source "$(dirname "${BASH_SOURCE[0]}")/build-wheel.bash"
			# TODO move building the wheels wo docker to RUN_DIRECT=1 build-wheel.bash ...

			# https://stackoverflow.com/questions/19331497/set-environment-variables-from-file-of-key-value-pairs
			# echo $(printf '%s' "$ENV" | awk -F= '{printf "%s=\"%s\" ",$1,$2 ;}' | xargs) bash -c "env"
			# cat <(printf '%s' "$ENV" | awk -F= '{printf "%s=\"%s\" ",$1,$2 ;}')
			VNE=$(printf '%s' "$ENV" | awk -F= '{printf "%s=\"%s\" ",$1,$2 ;}')
			# echo $VNE
			# env -S "$VNE" bash -c "env"
			(cd "${REPO_DIR}"; env -S "$VNE" bash -c "testing/linux/docker/build_on_docker.bash") #TODO refactor

			python3 -m pip uninstall -y eden-simulator

			WHEEL_TO_TEST=$(find "$BUILD_DIR/bin" -type f -name "eden_simulator-$VERSION-py3-none-*.whl")
			
			# TODO if dont run sphinx here...
			# "$(dirname "${BASH_SOURCE[0]}")/run-docs-with-wheel.bash" "$WHEEL_TO_TEST" TODO
			python3 -m pip install "$WHEEL_TO_TEST"
			python3 -m pip install $PIP_INSTALL_BUILD_DOCS_EXTRA -r "$REPO_DIR/docs/requirements.txt" # 
		fi
	fi
	# python3 -c "import eden_simulator; print(dir(eden_simulator));"
	
	if [ -z "$BUILD_IN_PLACE" ]; then
		rm -rf "$BUILD_DIR/docs" "$ARTIFACTS_DIR/*"
		# cp -r "$REPO_DIR/docs" "$BUILD_DIR/docs"
		cp -r "$REPO_DIR/." "$BUILD_DIR"
	fi
	
	# now build the docs!
	if [ -z "$DONT_RUN_SPHINX" ]; then # TODO a less awkward flag for readthedocs...
		python3 -m sphinx -T -E -W --keep-going -b html -d _build/doctrees -D language=en "${BUILD_DIR}/docs" $ARTIFACTS_DIR/html
	fi

else
	echo "Running \"reproducible\" docs on Docker, set RUN_DIRECT env var to override"
	source "$(dirname "${BASH_SOURCE[0]}")/docker/setup-options.bash"

	make -f "${REPO_DIR}/testing/linux/docker/Makefile" docker_build_docs_env

	# preferably use --env-file, otherwise: https://unix.stackexchange.com/questions/546053/convert-env-file-into-params-for-docker
	# or just run the same script inside docker

	bash "$(dirname "${BASH_SOURCE[0]}")/docker/sudo_docker.bash" \
		run -it --rm --mount type=bind,source=$(realpath ${REPO_DIR}),destination=/repo \
		--mount "type=bind,source=$(realpath $BUILD_DIR),destination=/build"  --user $DOCKER_USER_OR_ROOT \
		--mount "type=bind,source=$(realpath $ARTIFACTS_DIR),destination=/artifacts" \
		-e DEEP_IN_THOUGHT=1 -e RUN_DIRECT=1 -e BUILD_DIR=/build -e ARTIFACTS_DIR=/artifacts --workdir /build \
		eden-build-docs bash -c "set -e; bash /repo/testing/linux/build-docs.bash"

fi

# make this output only on the top level...
if [ -z "$DEEP_IN_THOUGHT" ]; then
	echo "Docs ready on $ARTIFACTS_DIR"
fi
