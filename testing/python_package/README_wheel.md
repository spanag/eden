# EDEN Python bindings + wheels

EDEN is a high-performance [NeuroML]( http://neuroml.org )-based neural simulator written in C++.
This wheel package contains Python bindings to use the simulator, as well as `win32`, `win_amd64`, `macosx_x86_64` and `manylinux1` standalone builds of the simulator.

## Installing

EDEN can be installed with one `pip` command:
```
pip install eden-simulator
```

The only non-`pip` requirement is for the GNU C compiler, Intel C compiler or Clang to be installed and available on `PATH` at run time.

## Usage

EDEN directly runs NeuroML v2 models of neural networks. 

The `eden_simulator` package exposes the `runEden` method, which takes the LEMS simulation file of the NeuroML model to be run as parameter:
```python
import eden_simulator
results = eden_simulator.runEden('<LEMS simulation file>.xml')
```
This interface returns the recorded trajectories specified in the simulation files in a Python dictionary, same as [`pyNeuroML`](https://pypi.org/project/pyNeuroML/) does with other simulation backends.

Thread-level parallelism can also be controlled with the `threads` keyword argument.

If other command-line arguments should be added, they can be specified as a list of strings in the optional parameter `extra_cmdline_args`.  
In case a specific instance of the EDEN executable is preferred, it can be selected with the optional parameter `executable_path`.  
Both options can be combined into a fully custom command to be run from Python: it can be passed as a list of strings in the optional parameter `full_cmdline`.

The recorded trajectories are also available relative to the working directory, as specified in the LEMS file.
Some `.gen.c` and `.gen.so` temporary files may also be generated in the working directory, these can be removed after the simulation is run.

For more information about the NeuroML model format, refer to http://neuroml.org .

The program's source code is available on GitLab: https://gitlab.com/neurocomputing-lab/Inferior_OliveEMC/eden .

This software is released as open-source, under the GNU General Public License version 3.

If you have any questions, suggestions or issues with running EDEN, please [submit a ticket]( https://gitlab.com/neurocomputing-lab/Inferior_OliveEMC/eden/-/issues ) on GitLab, or contact the developer at s.panagiotou@erasmusmc.nl .
