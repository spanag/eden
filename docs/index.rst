.. EDEN documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

================================
Welcome to EDEN's documentation!
================================

.. https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

**EDEN** (*E*\ xtensible *D*\ ynamics *E*\ ngine for *N*\ etworks) is a simulation program for `spiking neural networks <https://en.wikipedia.org/wiki/Spiking_neural_network>`_ that takes models described in `NeuroML <https://docs.neuroml.org/>`_ and generates the simulated behaviour of these networks.

..  TODO link TODO rate based?

It is best used as part of a neural modelling workflow, between the stage of *generating* the model to simulate, and *analysis* of the simulation's results.
	
To learn how to use EDEN, check out the `Quickstart <quickstart.ipynb>`_ with a classic Hodgkin-Huxley neuron, and the chapters of this guide for more about modelling and usage (coming soon).

.. ..., and the tutorials and our showcase of full-featured models with publication-ready figures.

Installing
**********

EDEN is most easily installed via :code:`pip`:

.. code::
	
	pip install eden-simulator

It can then be run from:
	- Python, as :code:`import eden_simulator as eden; eden.runEden('LEMS_<sim file>.xml')`;
	- the command line, as :code:`eden nml LEMS_<sim file>.xml`.

For more installation options, refer to the README `here <https://gitlab.com/c7859/neurocomputing-lab/Inferior_OliveEMC/eden/#installing>`_.

..
	.. code::
	

.. Note: this can't be added to the toctree as a reference to subsection on index.rst: https://github.com/sphinx-doc/sphinx/issues/2103 https://github.com/drdoctr/doctr/pull/285


.. Contact us 
.. .. include contact_us.rst

.. 
	code_cell
	usage
	
.. TODO use hidden toctree instead? who knows


	📽 Animations <anim>
	🎨 Render 3D <pyrender>

Table of contents
*****************

.. toctree::
	:maxdepth: 2
	:caption: Contents
	
	🚀️ Quickstart <quickstart>
	📖 Tutorial <tut_net>
	🐍️ Python API <python_api>


.. toctree::
	:caption: ⠀
	:maxdepth: 1
	
	🌟️ Examples <gallery>
	💌 Contact us <contact_us>	

.. toctree::
	:caption: More Examples
	:maxdepth: 1
	
	exa_lfp

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
