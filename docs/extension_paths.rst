================================
Extended LEMS paths
================================


A NeuroML model's variables and event sources can be designated in the form of `LEMS paths <https://docs.neuroml.org/Userdocs/Paths.html>`__.  These paths can be used for recording, among other uses.
However, there are some limitations of the paths currently used in NeuroML:

* LEMS is not enough to cover most aspects of spatial cells, and `discretisation <intro_spatial.ipynb>`__ is one of them.  Hence, there should be a way to refer to a specific location on a cell, in order to capture the individual compartments in all cases.
* The way that ``<Attachment>``\ s (that is, synapses and input sources) are accessed is ambiguous.  They are supposed to be numbered in the order that *mechanisms* were attached to the cells, *without distinction to the different component types, or logical groups of mechanisms attached* (for example, instances of the same synapse type that exist on the same cell are lumped together for all ``<projection>``\ s).  Even worse, the order that mechanisms are added is *unspecified*. 

Hence, EDEN also supports alternative forms to specify `precise locations <intro_spatial.ipynb#Specifying-a-location-on-a-cell>`__ on a spatially detailed cell, or a synapse or input mechanism, reliably.  These *extended LEMS paths* can be used to designate and record quantities and event sources alike, and EDEN's ``CustomSetup`` `extension <extension_customsetup.ipynb>`__ is also based on them.  They are written as follows:


LEMS paths for cell locations
*****************************
	
The typical way to access variables in a `spatially detailed <intro_spatial.ipynb>`__ cell is LEMS paths of the form:
	
* ``population/cell_id[/segment[.fractionAlong]]/``: (options continue from this path as follows)
	* ``v`` for membrane voltage
	* ``caConc`` for ``ca`` ion concentration
	* ``ca2Conc`` for ``ca2`` ion concentration (a separate ion species of the same atom)
	* ``biophysicalProperties/membraneProperties/``:
		* ``channel distribution/``:
		* ``iDensity`` or ``gdensity`` for the local current and conductance *per area* of the ion channel distribution
		* ``i`` or ``g`` for the current and conductance of the mechanism *over the whole pointed compartment* (recommended only in special cases)
		* ``channel mechanism/gate_name/q`` for gate variables
* Another way to specify a cell is ``population[cell_id]/...`` and the path follows the same way. 

The parts in square brackets are optional (if the part in enclosing brackets is used).  If ``segment`` is not specified, then it is assumed to be ``0`` (the soma, by convention).  If ``fractionAlong`` is not specified, then it is assumed to be ``0.5`` (the middle of the segment).

This way, any site on a neuron can be referred to with a LEMS path.  This is especially useful in conjunction with :py:func:`explain_cell <eden_simulator.experimental.explain_cell>` and  :py:func:`GetLemsLocatorsForCell <eden_simulator.experimental.GetLemsLocatorsForCell>`, and for fully recording over a cell, as seen in :doc:`intro_spatial` and other examples. 


LEMS paths for input list elements
**********************************

Instead of the :doc:`form <neuroml:Userdocs/Paths>` ``population/cell_id/synapses:mechanism name:serial number of attached instance/property``, which is supported by jLEMS (under certain assumptions), this is how quantities of ``input`` mechanisms can be accessed in EDEN as elements of their ``<inputList>``:

..
	LATER what about NEURON

.. code::
	
	input_list/0/variable_name
	input_list[0]/...


LEMS paths for synaptic projection elements
*******************************************

Instead of the LEMS form which may assign serial number in unpredictable ways, EDEN offers an alternative form similar to that for ``<inputList>``\ s.  After the syanpse instance and the mechanism's elements, a ``pre`` or ``post`` locator is added to select between either half of the synapse (for simple ``<projection>``\ s, only ``post`` is valid).  The pattern is then as follows:

.. code::
	
	projection/0/(pre or post)/variable_name
	projection[0]/(pre or post)/...

Using these forms, one can record or `otherwise access <extensions_pointers>`__ the variables of individual synapses or input sources, in the consistent order that they had in their respective ``<projection>``\ s or ``<inputList>``\ s.

LEMS paths for input streams
****************************

The same form can be further used to point to elements of `<𝙴𝚍𝚎𝚗𝚃𝚒𝚖𝚎𝚜𝚎𝚛𝚒𝚎𝚜𝚁𝚎𝚊𝚍𝚎𝚛> <extension_io.ipynb#Time-series-with-EdenTimeSeriesReader>`__\ s and `<𝙴𝚍𝚎𝚗𝙴𝚟𝚎𝚗𝚝𝚂𝚎𝚝𝚁𝚎𝚊𝚍𝚎𝚛> <extension_io.ipynb#Event-series-with-EdenEventSetReader>`__\ s and their properties, for `<𝚅𝚊𝚛𝚒𝚊𝚋𝚕𝚎𝚁𝚎𝚏𝚎𝚛𝚎𝚗𝚌𝚎> <extension_pointers.ipynb>`__\ s to point at or for recording:

* ``time series/instance[/column]`` (``column`` can be ignored if elements have just one)
* ``event  set/instance[/port]`` (``port`` can be ignored if elements have just one)
* And the ``group[instance]/...`` form can also be used as usual.
