# EDEN

_Extensible Dynamics Engine for Networks_

Parallel simulation engine for ODE-based network models.
Presently, implemented for NeuroML-based simulations of neural networks.

This software is released as open-source, under the GPL v3 license. Refer to LICENSE.txt file for more details.


## Quickstart
A Jupyter environment with EDEN and associated tooling pre-installed is available at https://github.com/spanag/eden-sim-jupyter-demo . 
Just click the Binder link on that page and you can use EDEN as shown on the bundled Python notebook.


## Installing

EDEN is currently available for the Linux platform.
The program can be installed through PyPI:
```sh
pip install eden-simulator
```
, or built from source with the provided Makefile. (The latter option is recommended for advanced uses, like MPI builds or system-wide installation)

The following software packages are required to build the source code:
- `gcc` compiler, or alternatively the `icc` compiler. Specifically, a compiler version that supports C++14.
- `flex`
- `bison` version 3.0 or later

Use environment variable `BUILD=release` or `BUILD=debug` to run a production or debugging build of the program executable, respectively. The executable will be available on `bin/eden.<build>.<compiler>.cpu.x`

If MPI is also installed, a MPI-enabled version of EDEN (with hybrid MPI/OpenMP parallelization) can be built, by running `make` with the `USE_MPI` flag.

Alternatively, Docker images with EDEN and an assortment of tools are available and can also be built, for containerized environments. The Dockerfiles are available on the `testing/docker` folder, and they can be built in the proper order through the Makefile in the folder.

If Docker is installed and accessible to the user, an automated testing suite can also be run to verify EDEN's results against the NEURON simulator's for deterministic models. Run `make test` to run the automated tests.


## Usage
EDEN directly runs NeuroML models of neural networks.  
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

Alternatively to the command line, the simulator can also be run within a Python program, if the `eden_simulator` Python package is installed.
The Python lines to run EDEN are then:
```python
import eden_simulator
results = eden_simulator.runEden('<LEMS simulation file>.xml');
```

This interface returns the recorded trajectories specified in the simulation files in a Python dictionary, same as pyNeuroML does with other simulation backends.
Thread-level parallelism can also be controlled with the `threads` optional argument.

For more information about the NeuroML model format, refer to http://neuroml.org .


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

