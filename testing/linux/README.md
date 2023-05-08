# Build scripts for Linux
The following instructions assume that the working directory contains the code repo under `eden/` subfolder.
This option can be changed by setting the environment variable `REPO_DIR`.


### Downloading the build tools

Make sure there are a few gigabytes of storage available, for the build tools.

In order to build EDEN on a typical Linux setup, certain tools have to be in place. Unfortunately, some tools must be installed on the system (thus modifying it) in order to build EDEN from source. These tools are:
- [Docker]( https://docs.docker.com/engine/install/ )
**OR** (coming soon):
- (Upcoming) [Homebrew]( http://brew.sh ) software package manager

The rest of the tools (namely: GNU Compiler Collection, recent flex and bison, Python3 and its wheel tooling) are installed as a Docker container.

A convenient script to set up all requirements in Docker is included in the same script used for building the release wheels and smoke tests using Docker, `eden/testing/mac/download-setup-requirements.bash`.

### Building EDEN

After the build tools are installed, EDEN can be built through the Makefile as usual, provided that the relevant tools are on `$PATH`. The `eden/testing/mac/setpath.bash` and scripts create a temporary directory with the necessary symlinks and set up `$PATH` to that purpose. The command line statements are thus like:
```
source eden/testing/mac/setpath.bash
cd eden
make eden <build options ...>
```

The command line to build a portable binary wheel is:
```sh
BUILD=release CFLAGS_extra="-static" WHEEL_PLAT_NAME_FILENAME=linux_x86_64 WHEEL_TARGET_PLAT=manylinux1_x86_64 make -j2 wheel
```
(The `manylinux1_x86_64` platform and `-j2` make parallelism are just suggestions.)
Be sure to test it on the corresponding manylinux container before publishing.

Finally, the Dockerized stages of building a EDEN wheel for the current CPU arch, and running smoke tests on that wheel, can be run via the script `eden/testing/mac/build-test-all-in-one.bash`.
After building and testing, the wheels are copied to a subdirectory of the working directory, `artifacts/`. This location can be overridden via the ARTIFACTS_DIR environment variable.

### Testing 

To run the full battery of tests requires validating against a 'gold' reference of what the results of the simulations should be; that reference is provided by the NEURON simulator.
Setting NEURON up automatically and without affecting the setup is done through a Docker on Linux environment.

Just run `make test` on the top level Makefile to run the full battery of tests.
