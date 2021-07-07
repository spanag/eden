# EDEN Python bindings + wheels

EDEN is a high-performance [NeuroML]( http://neuroml.org )-based neural simulator written in C++.
This wheel package contains Python bindings to use the simulator, as well as a `manylinux1`-compatible standalone build of the simulator.

## Installing

EDEN can be installed with one `pip` command:
```
pip install eden-simulator
```

The only non-`pip` requirement is for the GNU C compiler to be installed and available at run time.

## Usage

EDEN directly runs NeuroML v2 models of neural networks. 

The `eden_simulator` package exposes the `runEden` method, which takes the LEMS simulation file of the NeuroML model to be run as parameter:
```python
import eden_tools
results = eden_tools.runEden('<LEMS simulation file>.xml');
```
This interface returns the recorded trajectories specified in the simulation files in a Python dictionary, same as [`pyNeuroML`](https://pypi.org/project/pyNeuroML/) does with other simulation backends.

Thread-level parallelism can also be controlled with the `threads` keyword argument.

The recorded trajectories are also available relative to the working directory, as specified in the LEMS file.
Some `.gen.c` and `.gen.so` temporary files may also be generated in the working directory, these can be removed after the simulation is run.

For more information about the NeuroML model format, refer to http://neuroml.org .

The program's source code is available on GitLab: https://gitlab.com/neurocomputing-lab/Inferior_OliveEMC/eden .

This software is released as open-source, under the GNU General Public License version 3.

If you have any questions, suggestions or issues with running EDEN, please [submit a ticket]( https://gitlab.com/neurocomputing-lab/Inferior_OliveEMC/eden/-/issues ) on GitLab, or contact the developer at s.panagiotou@erasmusmc.nl .
