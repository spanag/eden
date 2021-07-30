# Build scripts for Windows

The following instructions assume that the working directory contains the code repo under `eden/` subfolder, and the build tools will be installed in a new sibling folder, `tooling/`. These options can be changed by setting the environment variables `REPO_DIR` and `TOOLING_DIR` respectively.


### Downloading the build tools

Make sure there are around 2 gigabytes of storage available, for the build tools; even more if running the tests on Docker.

In order to build EDEN on a typical Windows setup, certain tools have to be in place. Fortunately, all these tools can be just placed on a directory, without the need to alter the Windows setup's list of installed programs, PATH or Registry entries, or other persistent changes.
A convenient script to set up all requirements is `eden\testing\windows\download-setup-requirements.bat`.

This downloads and extracts all requirements into a child directory `%TOOLING_DIR%` of the working directory, thus the script is preferably run on a working directory outside the source tree.

To begin downloading the required tools, `wget` is downloaded using the standard PowerShell facility for downloading files; if PowerShell is not present, the relevant line script should be changed to use an existing tool, or it could be commented out with `wget` being provided manually by the user.

This script is provided as a convenience only, since the URL's it references are valid only for the time the script was last edited. It would be best if the binaries for the tools were downloaded once, verified, and then pulled from a secure repository; feel free to modify the script to that end.


### Building EDEN

After the build tools are in place, EDEN can be built through the Makefile as usual, provided that the tools are on `%PATH%`. The `eden\testing\windows\setpath-i686.bat` and `eden\testing\windows\setpath-amd64.bat` scripts set up `%PATH%` in order to target 32-bit x86, or x86-64 respectively. The command line statements are thus like:
```cmd
eden\testing\windows\setpath-i686.bat
cd eden
make eden <build options ...>
```

The provided scripts `eden\testing\windows\build-eden-i686.bat` and `eden\testing\windows\build-eden-amd64.bat` build EDEN and package the Python wheels for the respective CPU architectures.


### Testing 

To run the full battery of tests requires validating against a 'gold' reference of what the results of the simulations should be; that reference is provided by the NEURON simulator.
Setting NEURON up automatically and without affecting the setup is a work in progress; for now the full tests are run on a Docker on Linux environment. A Docker environment for EDEN on Windows is also in development.

Instead, for the moment, a simple smoke test is run on Windows. Code paths are largely the same among platforms, the only practical difference is in generating code and invoking a compiler at runtime, which is what is being tested.

The script to run the tests that are run in the moment for Windows, is `eden\testing\windows\run-basic-tests.bat`. Python has to be in `PATH`, with `eden-simulator` installed.  
The script to run the tests for a specific wheel is `eden\testing\windows\run-tests-on-wheel.bat`, with a command-line argument of which wheel to test. Python has to be in `PATH`. 

Finally, all the above stages of setting up tooling, building, and basic tests can be run on a fresh Windows system via the `eden\testing\windows\setup-build-test-all-in-one.bat`.

The wheels are copied to a subdirectory of the working directory, `artifacts/`. This location can be overridden via the ARTIFACTS_DIR environment variable.
