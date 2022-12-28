#!/bin/bash
set -e

# get brew path from here
source "$(dirname "${BASH_SOURCE[0]}")/setup-options.bash"
PATH="$HOMEBREW_BIN:$PATH"

if [[ -z $HOMEBREW_INSTALL_LOCAL ]]; then
    # if brew is not found:
    if [ "$(which brew)" == "" ]; then
        # Set up Homebrew, on the default path
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        brew analytics off
    fi
else
    # https://docs.brew.sh/Installation "Untar anywhere"
    mkdir homebrew && curl -L https://github.com/Homebrew/brew/tarball/master | tar xz --strip 1 -C $HOMEBREW_INSTALL_LOCAL_DIR
    eval "$($HOMEBREW_INSTALL_LOCAL_DIR/bin/brew shellenv)"
    brew analytics off
    brew update --force --quiet
fi

# flex is not tagged, apparently?
EXTRA_BREW_PACKAGES=
if  [ "$BUILD_TARGET_ARCH" = "arm64" ]; then
    EXTRA_BREW_PACKAGES=" hdf5@1.12"
fi
brew install gcc@11 flex bison@3.8 python@3.9 $EXTRA_BREW_PACKAGES

TEMP_LINK_DIR=$(pwd)/build_aliases
python3 -m venv $TEMP_LINK_DIR/venv
source $TEMP_LINK_DIR/venv/bin/activate

python3 -m pip install -U pip setuptools wheel delocate twine

echo "Requirements installed"
