#!/bin/bash
set -e

# Set up Homebrew, on the default path
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
# brew analytics off

# Or alternatively, install on any old path on, say, home. Though then, gcc will have to be built from source, which takes a while !
# mkdir homebrew && curl -L https://github.com/Homebrew/brew/tarball/master | tar xz --strip 1 -C homebrew
# PATH=$(pwd)/homebrew/bin:$PATH

brew install gcc@11 flex bison python3

python3 -m pip install -U pip setuptools auditwheel twine

