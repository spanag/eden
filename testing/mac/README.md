# Build scripts for Mac OS X

The following instructions assume that the working directory contains the code repo under `eden/` subfolder.
This option can be changed by setting the environment variable `REPO_DIR`.


### Downloading the build tools

Make sure there are a few gigabytes of storage available, for the build tools.

In order to build EDEN on a typical Mac setup, certain tools have to be in place. Unfortunately, some tools must be installed on the system (thus modifying it) in order to build EDEN from source. These tools are:
1. [Command Line Developer Tools for Mac]( https://developer.apple.com/downloads/index.action?=command%20line%20tools )
2. [Homebrew]( http://brew.sh ) software package manager

The rest of the tools (namely: GNU Compiler Collection,  recent flex and bison, Python3) are installed as Homebrew packages.  

A convenient script to set up all requirements is `eden/testing/mac/download-setup-requirements.bash`.
This installs Homebrew on its default (system-wide) location if it is not already available, gets the command line developer tools also installed in the process, and gets the necessary Homebrew packages installed.

This script is provided as a convenience only; feel free to modify the script to your specific needs (such as when the necessary tooling is already installed, by a different system).


### Building EDEN

After the build tools are in place and gcc and g++ available on PATH are the real deal (instead of Apple-built `clang`), EDEN can be built through the Makefile as usual, provided that the relevant tools are on `$PATH`. The `eden/testing/mac/setpath.bash` and scripts create a temporary directory with the necessary symlinks and set up `$PATH` to that purpose. The command line statements are thus like:
```
source eden/testing/mac/setpath.bash
cd eden
make eden <build options ...>
```

The provided script `eden/testing/mac/build-wheel.bash` builds EDEN and packages the `eden-simulator` Python wheel for the machine's CPU architecture (or the simulated one, if it is run on Rosetta).


### Testing 

To run the full battery of tests requires validating against a 'gold' reference of what the results of the simulations should be; that reference is provided by the NEURON simulator.
Setting NEURON up automatically and without affecting the setup is a work in progress; for now the full tests are run on a Docker on Linux environment. A Docker environment for EDEN on macOS is also in development.

Instead, for the moment, a simple smoke test is run on macOS. Code paths are largely the same among platforms, the only practical difference is in generating code and invoking a compiler at runtime, which is what is being tested.

The script to run the tests that are run in the moment for macOS, is `eden/testing/mac/run-basic-tests.bash`. Python has to be in `PATH`, with `eden-simulator` installed.  
The script to run the tests for a specific wheel is `eden/testing/mac/run-tests-on-wheel.bash`, with a command-line argument of which wheel to test. Python has to be in `PATH`. 

Finally, the stages of building a EDEN wheel for the current CPU arch, and running basic tests on that wheel, can be run via the script `eden/testing/mac/build-test-all-in-one.bash`.
After building and testing, the wheels are copied to a subdirectory of the working directory, `artifacts/`. This location can be overridden via the ARTIFACTS_DIR environment variable.
