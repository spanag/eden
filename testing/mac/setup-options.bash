#!/bin/bash
# Some configuration options are placed here to be accessible from both setup-install-requirements and also regular build-test scripts,
# basically because Homebrew takes too long to install when placed in an uncommon prefix, and also it is not on PATH in the new M1 Mac's.

# LATER make this a supplied parameter if multiarch is needed
BUILD_TARGET_ARCH=$(uname -m)
# LATER check if it's an identified arch. for now assume either x86_64 or arm64

# Alternatively, install on any old path on, say, cwd (above eden repo). Though then, gcc will have to be built from source, which takes a while !
# HOMEBREW_INSTALL_LOCAL=1
# HOMEBREW_INSTALL_LOCAL_PREFIX=$(pwd)

if [[ -z $HOMEBREW_INSTALL_LOCAL ]]; then
    # First, try finding Homebrew just in case it is on PATH for whatever reason
    if [ "$(which brew)" != "" ]; then
        HOMEBREW_BIN=$(dirname $(which brew))
    else
        # global install it is, then, use the Homebrew default
        # but that is different per target arch ...
        if [ $BUILD_TARGET_ARCH == "arm64" ]; then
            HOMEBREW_BIN=/opt/homebrew/bin
        else
            HOMEBREW_BIN=/usr/local/bin
        fi
    fi
    # echo $HOMEBREW_BIN
else
    # let's hardcode the "homebrew" dir name for now since it's probably for a self contained build env
    HOMEBREW_INSTALL_LOCAL_DIR=$HOMEBREW_INSTALL_LOCAL_PREFIX/homebrew
    HOMEBREW_BIN=$HOMEBREW_INSTALL_LOCAL_DIR/bin
fi
