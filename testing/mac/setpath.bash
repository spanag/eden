#!/bin/bash

# Get brew path from here
source "$(dirname "${BASH_SOURCE[0]}")/setup-options.bash"
PATH="$HOMEBREW_BIN:$PATH"

# Create a temporary folder, relative to cwd, to store links to the binaries to be used
# something more reasonable LATER

TEMP_LINK_DIR=$(pwd)/build_aliases

mkdir -p $TEMP_LINK_DIR

# These are based on Homebrew with a specific version, for now
# TODO select version of gcc
HOMEBREW_PATH=$(brew --prefix)
HOMEBREW_BIN="$HOMEBREW_PATH/bin"
ln -fs "$HOMEBREW_BIN/gcc-11" "$TEMP_LINK_DIR/cc"
ln -fs "$HOMEBREW_BIN/gcc-11" "$TEMP_LINK_DIR/gcc"
ln -fs "$HOMEBREW_BIN/g++-11" "$TEMP_LINK_DIR/g++"
ln -fs "$HOMEBREW_BIN/c++-11" "$TEMP_LINK_DIR/c++"
ln -fs "$HOMEBREW_BIN/cpp-11" "$TEMP_LINK_DIR/cpp"

ln -fs "$HOMEBREW_PATH/opt/bison/bin/bison" "$TEMP_LINK_DIR/bison"
ln -fs "$HOMEBREW_PATH/opt/flex/bin/flex"   "$TEMP_LINK_DIR/flex"

PATH="$TEMP_LINK_DIR:$PATH"

source $TEMP_LINK_DIR/venv/bin/activate

# may also need HDF5 install present on M1 Mac's
# set it here because HOMEBREW_PATH is used
if  [ "$BUILD_TARGET_ARCH" = "arm64" ]; then
    export HDF5_DIR=$HOMEBREW_PATH
fi
