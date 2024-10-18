.. EDEN user manual top level file

================================
User's Guide
================================

This chapter introduces how to use EDEN, its modelling language and simulation capabilities, and it shows how each feature can be used to build many types of neural models.

.. can't put comments on toc links, thus we'll write this note here https://stackoverflow.com/questions/38836458/annotated-sphinx-toctree

:ref:`/quickstart.ipynb` shows how to run EDEN with ready-made files; to learn modelling in NeuroML, start from :ref:`/intro_neuroml.ipynb`.


User guide contents
-------------------

.. toctree::
	:maxdepth: 1
	:caption: Introducing NeuroML
	
	Quickstart <quickstart>
	NeuroML basics <intro_neuroml>
	Detailed cells with NeuroML <intro_spatial>
	LEMS basics <intro_lems>
	
.. NEXT replace extensions-list with intro

:hidden:`Extensions list`
-------------------------

.. toctree::
	:maxdepth: 1
	:caption: Beyond NeuroML: EDEN-specific extensions
	
	extension_paths
	extension_customsetup
	example_spatial_customsetup
	extension_pointers
	extension_io
	extension_multiflux
	extension_writable

:hidden:`Examples`
------------------

.. toctree::
	:caption: Model Examples
	:maxdepth: 1
	
	Putting it all together: A network of detailed cells <tut_net>
	exa_lfp
	example_pong

..
	example_imposed_field
	example_robot
	example_mea
	example_matlab

..  toctree
	:maxdepth: 2
	:caption: Reference

.. NEXT html?

See also:

* 🐍️ :doc:`python_api`
* 🌟️ :doc:`Examples gallery <gallery>`

.. 	
	📽 Animations <anim>
	🎨 Render 3D <pyrender>

