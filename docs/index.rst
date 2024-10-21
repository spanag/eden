.. EDEN documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

================================
Welcome to EDEN's documentation!
================================

.. https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

.. include:: foreword.rst

.. include:: installing.rst

..
	.. code::
	

.. Note: this can't be added to the toctree as a reference to subsection on index.rst: https://github.com/sphinx-doc/sphinx/issues/2103 https://github.com/drdoctr/doctr/pull/285

.. NEXT use hidden toctree instead? who knows


	📽 Animations <anim>
	🎨 Render 3D <pyrender>

.. rubric:: Table of contents
	:name: indextoc
	:heading-level: 2

.. toctree::
	:maxdepth: 2
	:caption: Contents
	:titlesonly:
	
	🚀️ Quickstart <quickstart>
	📖 User's Guide <user_guide>
	🐍️ Python API <python_api>

.. toctree::
	:caption: ⠀
	:maxdepth: 1
	
	🌟️ Examples <gallery>
	 ❓ FAQ <faq>
	💌 Contact us <contact_us>	

..
	Intro
		Installing, Contact us, Credits
	
	Part A - User's Guide
		Intro to NeuroML
		More than NeuroML
		Modelling examples
		Reference? thips and tricks?
			rendering, animation?
	
	Part B: Hacker's guide
		Theory
		Maintenance?
		Examples
	
	appendix - Reference
		FAQ
		Python API
		C++ API?

.. rubric:: Indices and tables
	:heading-level: 3

* :ref:`genindex`
* :ref:`search`

.. 
	* :ref:`modindex`
