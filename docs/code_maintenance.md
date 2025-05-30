# Software Maintenance

This section concerns the standard EDEN distribution's codebase and build process. It is useful for both maintainers of the codebase, as well as who wish to interface with it. 

<!-- Building and testing -->

## Testing process

For Linux: with Docker available, run `make test`.
Currently, proper regression testing is limited to Linux for want of time.  A special Docker image is built and run, to bundle all necessary software and isolate the test environment from the user's own.
Testing per platform can be added, if the other related software (NEURON, NeuroML tools, build tools) and their dependencies can also run. (binaries are often harder to find for non-PC CPUs, and such)
A basic smoke test is done at the end of the build process, for all platforms; this is not enough but it can catch many platform specific errors. (Isolation between the build and the smoke-test environment is a work in progress for OSX, Linux, Windows ordered by severity.)

The combinations of interacting mechanisms are way too many to completely cover with tests. To keep bugs in check, make sure to maintain regular, clean-cut interfaces between code modules. Minimise code duplication, re-use elements and use strong typing to increase the likelihood of catching bugs, and the effectiveness of bug fixes. (A bit of duplicated code structure may often be better than a single routine pockmarked with conditionals; exercise judgement.)

Keep modules tight; using strictly one source file per class may be counter-productive. (See also: Java and the usefulness of inner classes.)

## Build process

The primary incremental build tool should be GNU `make`.  Although newer build tools like CMake and Meson are more elegant and automagical, this often backfires on custom exotic systems like e.g. supercomputers.  In such cases, `make`'s direct and predictable way of issuing commands is much easier to control and adapt. Meanwhile, the older tools like autoconf/automake are not so important given the environment convergence going on and the limited complexity of this program, and besides they may also behave in magical ways.

A good example for how a Makefile should be is the [Makefile](https://github.com/git/git/blob/master/Makefile) of the Git program.

This makefile is primarily designed for the GNU Compiler Collection, as it's good for all common platforms.  Thanks to user-facing convergence, the compiler options are also largely valid for the clang compiler as well, and likely the new, Clang based Intel compiler.

If it is desired to build with Visual Studio, GNU Make can still be used (look it up if interested). The compiler flags and linkage will need a full overhaul for VC++ though, and edits to the C++ source may also be in order (We try to keep the code standards compliant, but each compiler has multiple conflicting opinions on what is valid C++).


### Build environment

- On Linux, builds are based on `manylinux1` for regular PCs, and `manylinux2014` for most other CPU architectures. A version bump may be needed for newer C++ versions and facilitated by the update treadmill in the userbase.
For Windows and OSX, the required set of build/test tools is assembled using the `testing/<os>/download-requirements` scripts.
- On Windows, the software binaries are gathered from various sources. On OSX, they are downloaded and installed with Homebrew; check the install script to separate the home environment from Homebrew's default. (Though all binaries will have to be rebuilt in that case; consider a virtual machine for easy isolation.)
- In principle, automated testing can be set up easily.  It will just take some time to set up the tools for testing, the same way it's been done for building.

### Release binaries
Binaries are distributed in the Python wheel format, because `pip install eden-simulator` is more convenient and less error prone than manually selecting the correct arch for the target machine and unpacking.  [Python bindings]( python_api.rst ) are of course included with the wheels.

Wheels are generated using scripts found at `testing/<OS>/build-test-all-in-one` for each operating system (where they are supposed to be built on). **NB**: only basic "smoke tests" are being run, full testing is done on Linux with `make test` [for now]( #testing-process ).

Linker issues for release binaries are mitigated by the usual DLL bundling tools, `delocate` and `auditwheel` for Mac and Linux respectively.  `delvewheel` can also be used for Windows, but these builds can remain static for now.  (Be aware that careless, automatic DLL bundling can violate code licensing; watch out for unintended dependencies.)


<!-- Makefile TODO
It may be augmented with scripts in a language than is popular and easy to obtain (eg Bash/CMD, Python)
options... -->


<!-- Git management TODO -->

<!-- Release process
one test run
TBD TODO -->

<!-- Docs build process TODO -->


:::{spelling:word-list}
autoconf
automake
:::
