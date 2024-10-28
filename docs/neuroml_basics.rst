NeuroML basics
______________

NeuroML is (as of this writing) the most extensive data format to describe models of spiking neural networks (SNN's) (and generally dynamic neural networks), outside the context of a specific simulation system.  NeuroML version 2 (i.e. the latest one) relies on the LEMS language for defining new (ODE and event-based) dynamics for the various parts that constitute a neural network.  EDEN runs models thet are described in this NeuroML-v.2 (plus LEMS) data format.

This part introduces the basic concepts of NeuroML (in its second version) and LEMS, and how these can be applied to construct full-featured models of dynamic neural networks.

.. toctree::
	:maxdepth: 2
	
	Quickstart <quickstart>
	NeuroML basics <intro_neuroml>
	Detailed cells with NeuroML <intro_spatial>
	LEMS basics <intro_lems>
