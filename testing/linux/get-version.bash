#!/bin/bash

source "$(dirname "${BASH_SOURCE[0]}")/default-repo-path.bash"

if [ ! -f "$REPO_DIR/VERSION" ]; then
	echo "$REPO_DIR/VERSION" not found, "$REPO_DIR" is probably not the root of the EDEN repo !
	exit 2
fi
VERSION=$(cat "$REPO_DIR/VERSION")
# echo $VERSION
