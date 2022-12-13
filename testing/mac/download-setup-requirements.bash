#!/bin/bash
set -e

# Set up Homebrew, on the default path
/bin/bash -c "brew --version || $(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew analytics off

# Or alternatively, install on any old path on, say, home. Though then, gcc will have to be built from source, which takes a while !
# mkdir homebrew && curl -L https://github.com/Homebrew/brew/tarball/master | tar xz --strip 1 -C homebrew
# PATH=$(pwd)/homebrew/bin:$PATH

brew install gcc@11 flex bison python3

TEMP_LINK_DIR=$(pwd)/build_aliases
python3 -m venv $TEMP_LINK_DIR/venv
source $TEMP_LINK_DIR/venv/bin/activate

python3 -m pip install -U pip setuptools wheel delocate twine

echo "Requirements installed"
