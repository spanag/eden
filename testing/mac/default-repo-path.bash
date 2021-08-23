#!/bin/bash


if [ -z "$REPO_DIR" ]; then
	REPO_DIR=$(builtin cd $(dirname "${BASH_SOURCE[0]}")/../../; pwd)
	echo REPO_DIR was not defined, set to $REPO_DIR
fi
