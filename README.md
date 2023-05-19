# EDEN

_Extensible Dynamics Engine for Networks_

Parallel simulation engine for ODE-based network models.
Presently, implemented for NeuroML-based simulations of neural networks.

This software is released as open-source, under the GPL v3 licence. Refer to LICENCE.txt file for more details.


## Quickstart

The easiest way to get started is to get the Python package: `pip install eden-simulator`
and then invoke EDEN from Python:
```python
from eden_simulator import runEden
sim_results = runEden('<LEMS simulation file>.xml') # replace filename with your own
```

For a demo run of NeuroML networks and performance comparison with NEURON:
A Jupyter environment with EDEN and associated tooling pre-installed is available at https://github.com/spanag/eden-sim-jupyter-demo . 
Just click the Binder link on that page and you can use EDEN as shown on the bundled Python notebook.


## Installing

EDEN is currently available for the Microsoft Windows™, Apple macOS™ and Linux operating systems.
The program can be installed through PyPI:
```sh
pip install eden-simulator
```
, or built from source with the provided Makefile. (The latter option is recommended for advanced uses, like MPI builds, exotic platforms, or re-packaging)

### A compiler must be available at run time
Since EDEN generates code tailored to the model to be run each time, a C compiler must be available and accessible from `PATH`. The GNU and Intel compilers are currently supported for this purpose.
If the compiler to be used (`gcc` by default) is not accessible on `PATH`, the simulation fails at setup and an error message with suggested actions is generated.

On Windows setups, there is usually no compiler installed. A tested version of GCC can be downloaded for use with EDEN, from:  
https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/8.1.0/threads-posix/sjlj/i686-8.1.0-release-posix-sjlj-rt_v6-rev0.7z (for 32-bit EDEN)  
https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/8.1.0/threads-posix/seh/x86_64-8.1.0-release-posix-seh-rt_v6-rev0.7z (for 64-bit EDEN)  
Unpack the file anywhere, and add the unpacked `<path ...>\bin` directory to EDEN's PATH. (More on how to do this [below](#setting-path-to-include-a-compiler))

On macOS setups, a compatible compiler can be installed with the [Command Line Developer Tools for Mac]( https://developer.apple.com/downloads/index.action?=command%20line%20tool ). Follow the link for instructions, or simply run `xcode-select --install` from the macOS terminal and follow the GUI prompts.

On Linux setups, GCC is commonly installed; if it is not, refer to your distribution's documentation on how to install it (or to respective vendors for other compilers).

In all cases, make sure that the compiler targets the same architecture as the executable, see also [Usage](#usage) for more details.


## Usage

### From the command line
EDEN directly runs NeuroML models of neural networks. (For more information about the NeuroML model format, refer to http://neuroml.org ) 
It can be run from the command line, with the following command referencing the LEMS simulation file of the NeuroML model to be run:  
```
<executable> nml <LEMS simulation file> 
```  
For example,  
``` 
bin/eden.release.gcc.cpu.x nml examples/LEMS_NML2_Ex25_MultiComp.xml
``` 

The recorded trajectories will be available relative to the working directory, as specified in the LEMS file.
Some `.gen.c` and `.gen.so` temporary files may also be generated, these can be removed after the simulation is run.

Per-thread parallelism can be adjusted through the `OMP_NUM_THREADS` environment variable.

### From Python
Alternatively to the command line, the simulator can also be run within a Python program, if the `eden_simulator` Python package is installed.
The Python lines to run EDEN are then:
```python
import eden_simulator
results = eden_simulator.runEden('<LEMS simulation file>.xml') # replace filename with your own
```
This interface returns the recorded trajectories specified in the simulation files in a Python dictionary, same as pyNeuroML does with other simulation backends.

Thread-level parallelism can also be controlled with the `threads` optional argument.

If other command-line arguments should be added, they can be specified as a list of strings in the optional parameter `extra_cmdline_args`.  
In case a specific instance of the EDEN executable is preferred, it can be selected with the optional parameter `executable_path`.  
Both options can be combined into a fully custom command to be run from Python: it can be passed as a list of strings in the optional parameter `full_cmdline`.

### Setting `PATH` to include a compiler

When running EDEN, make sure that the selected compiler (`gcc` by default) is available on PATH and has the same bitness and architecture as the build of EDEN in use. This is a concern on Windows (where both 32 and 64 bit programs are common) and on certain HPC clusters where the login and job nodes may run different instruction sets.

The `PATH` can be set this way on the Windows command line:
```cmd
path "full_path_to_compiler_executable";%PATH%
eden ...
```
And this way, in Unix shells (Linux, OSX, AIX, etc...)
```sh
PATH="full_path_to_compiler_executable":$PATH eden ...
```
If using Python, `PATH` can be set as follows:
```python
os.environ["PATH"] = <path to compiler executable> + os.pathsep + os.environ["PATH"]
runEden(...)
```


## Building from source
The following software packages are required to build the source code:
- `gcc` compiler, or alternatively the `icc` compiler. Specifically, a compiler version that supports C++14.
- `flex` version 2.6 or later
- `bison` version 3.0 or later
- (optional) `python3` with `setuptools` and optionally `wheel` to build the Python wrapper package (see also the next section).

for building on Linux, shell scripts are available on the [testing/linux](testing/linux) directory.
The necessary tools for building are pre-installed or they can be easily installed; consult your distribution's reference on installing essential build tools. 
Alternatively, if `docker` is installed and available to the user, EDEN can be built without affecting the user setup by just running `bash eden/testing/linux/build-test-all-in-one.bash` outside the source tree. The executable and Python wheel can then be found in `eden/bin` directory.

For building on Windows, batch scripts are available on the [testing/windows](testing/windows) directory for installing all the necessary build tools and libraries, setting the shell's PATH to use them and building the standalone executable or Python wheel, all without affecting the system or user setup.

For building on macOS, shell scripts are available on the [testing/mac](testing/mac) directory for installing all the necessary build tools and libraries, setting the shell's PATH to use them and building the standalone executable or Python wheel.
*Note* The `testing/mac/download-setup-requirements.bash` script installs Command Line Developer Tools for Mac, Homebrew and various Homebrew packages in the system, in the process installing the required tooling.
*Note 2* The Apple developer tools are not suitable to build EDEN, because they lack OpenMP support which is required by EDEN. Instead, we recommend the developer tools installed by `download-setup-requirements.bash`.

Before attempting to build manually on Windows or macOS, source the `setpath` script for the platform you are building against so that the necessary tools are on PATH. (See platform-specific instructions above.)

When running `make eden` (or `make wheel`, see below) manually, use environment variable `BUILD=release` or `BUILD=debug` to run a production or debugging build of the program executable, respectively. The executable will be available on `bin/eden.<build>.<compiler>.cpu.x` (`.exe` on Windows)

### Building Python wheels
Beside the program itself, EDEN can also be built as a Python package, which offers more integrated interface to the program. The python package is created in the installable `.whl` (wheel) format. 

There are two options to build a Python wheel:
- The first is as a standalone wheel containing the EDEN executable (hence the wheel is specific to one OS and processor type).  This is the type of wheels available through `pip install`.
- The second option is as a 'hollow' wheel which works everywhere, but relies on EDEN already being available on PATH through different means. This type is useful in classes where a special (e.g. custom-built) version of EDEN is preferred to the generic version of 'standalone' wheels.

Both types of wheel can be built as the respective targets `wheel` and `hollow_wheel` of the Makefile. The resulting `.whl` files are located in the paths `testing/sandbox/{wheel, wheel_hollow}/dist`.

### Building for MPI
If MPI is also installed, a MPI-enabled version of EDEN (with hybrid MPI/OpenMP parallelization) can be built, by running `make` with the `USE_MPI` flag. This configuration has been tested with standard MPICH on Linux; consult your HPC cluster's documentation for specific details on the process for MPI-enabled builds.

## Docker images
Alternatively, Docker images with EDEN and an assortment of tools are available and can also be built, for containerized environments. The Dockerfiles are available on the `testing/linux/docker` folder, and they can be built in the proper order through the Makefile in the folder.
The Docker images are Linux-native but they can as well run on Windows and MacOS through [Docker Desktop]( https://www.docker.com/products/docker-desktop/ ).

If Docker is installed and accessible to the user, an automated testing suite can also be run to verify EDEN's results against the NEURON simulator's for deterministic models. Run `make test` to run the automated tests.


## Dockerfile

A demonstration-ready Docker image is available, using the Dockerfile on this top-level directory.

Beside EDEN, the demo image also contains the following:

- Jupyter Notebook
- NEURON simulator package, also for neural-network simulations
- pyNeuroML/jNeuroML, to easily run NeuroML-based simulations and access results through Python notebooks

To use, first build the Docker image: 

	docker build -t eden_demo .

To run the Docker image:

	docker run -it --name eden_demo_book -v ~/user_work_folder:/home/jovyan/work -p 8888:8888 --rm eden_demo

(`user_work_folder` can be replaced by any folder the user is working on; this is where user files will be read from, and work done using the image will be stored)

Then, use the URL displayed on the terminal to work with Jupyter and the included software.


## Acknowledgements

This research is supported by the European Commission Horizon2020 Framework Programme Projects EXA2PRO (Grant Agreement No. 801015) and EuroEXA (Grant Agreement No. 754337).

EDEN uses third-party open source libraries. The code of the included libraries, and further information, can be found in the `thirdparty/` top-level directory.


## Contact
For support and questions, contact Sotirios Panagiotou < s.panagiotou@erasmusmc.nl > and Dr. Christos Strydis < c.strydis@erasmusmc.nl > .

