'''
Additional facilities to use EDEN's experimental features from within Python.

**NOTE:** As the experimental features are still in design, parts of this API may change between versions!
'''

import json
import numpy as np
from .run_sim import _invoke_eden

def explain_cell( nml_file, *, verbose = False, threads=None,
	extra_cmdline_args=None, executable_path=None, full_cmdline=None,
	**kwargs
):
	'''
	Get more information about the structure and properties of physically-modelled cells.
	
	
	Parameters
	----------
	
	nml_file : str
		The main NeuroML file (and the ones it <include>s) to consider.
	
	verbose : bool, optional
		Add some more console output when running the simulator.
		
	threads : int, optional
		The number of threads to run processing with, or `None` for automatic selection. (Does not have an effect in this case yet.)
		
		⠀
		
		*The following parameters are less commonly used*:
	
	
	Returns
	-------

	info : dict[str, dict]
		The requested information, keyed by cell type name.
	
	
	Other Parameters
	----------------
	
	
	extra_cmdline_args : list[str], optional
		Additional command line arguments to pass to EDEN.  
		Refer to the EDEN user's manual, "Command line API" for more details.
		
	full_cmdline : list[str], optional
		Specify the exact argv to be passed to EDEN.  
		Refer to the EDEN user's manual, "Command line API" for more details.
	
	executable_path : str, optional
		Select a specific executable of EDEN to run the simulation with. Useful for custom installations and uses of EDEN.
	
	
	Notes
	-----
	
	*For physical (i.e. not artificial) cells:*
	
	Each keyed entry of `info` may have the following dict entries\: 
	
	``'comp_parent'``: ndarray[(n_comps), int]
		The tree parent of each comparent. The first value, that represents the `root` (as modelled) of the neuron, is -1.
	
	``'comp_start_pos'``: ndarray[(n_comps, 3), float]
		The position of the most ``proximal`` end of each compartment.
	
	``'comp_end_pos'``: ndarray[(n_comps, 3), float]
		The position of the most ``distal``   end of each compartment.
	
	``'comp_midpoint'``: ndarray[(n_comps, 3), float]
		The 3-D position coordinates of the midpoint of each compartment, relative to the morphology's origin.  
		*NB:* Due to curvature it may not lay exactly between `comp_start_pos` and 'comp_end_pos', or on the line joining them.
	
	``'comp_midpoint_segment'``: ndarray[(n_comps), int]
		The NeuroML <segment> containing the midpoint of each compartment. 
		
	``'comp_midpoint_fractionAlong'``: ndarray[(n_comps), float]
		The 0-to-1, proximal-to-distal, NeuroML <segment> position where each compartment lies on its `comp_midpoint_segment`.  
	
	``'comp_length'``: ndarray[(n_comps), float]
		The length of the neurite cable belonging to each compartment.  
		*NB:* Due to curvature and explicit discontinuities this may be more or less than the straight-line distance between `comp_start_pos` and `comp_end_pos`.
	
	``'comp_path_length_from_root'``: ndarray[(n_comps), float]
		The distance from the neuron's root tracing back the anatomical path, in microns.  Some biophysical attribute track this distance.
	
	``'comp_area'``: ndarray[(n_comps), float]
		The total membrane area per compartment, in `μm²`.
	
	``'comp_volume'``: ndarray[(n_comps), float]
		The total enclosed volume per compartment, in `μm³`.
	
	``'comp_capacitance'``: ndarray[(n_comps), float]
		The total membrane capacitance per compartment, in `pF`.
	
	``'comp_conductance_to_parent'``: ndarray[(n_comps), float]
		The electrical cytosolic conductance betwen each comparent and its tree parent, in `nS`. The first value for the tree root is zero (not applicable).
	
	``'segment_groups'``: dict[str, dict]
		Details about each NeuroML segment group, keyed by name.
		
		Each keyed entry of `segment_groups` may have the following dict entries\: 
		
		* ``'comps'``: The numbered compartments included in the segment group.
	
	``'mesh_vertices'``: ndarray[(n_verts, 3), float]
		The vertex coordinates of the 3-D mesh representing the neuron, in microns.
	
	``'mesh_faces'``: ndarray[(n_faces, 3), int]
		The triangles of *0-indexed vertices* forming the neuron's 3-D mesh.
		
	``'mesh_comp_per_face'``: ndarray[n_faces, int] (optional)
		The neuron compartment that each triangle of the mesh maps to.  
		
		Useful for selectively colouring relevant regions of the neuron, and displaying the neuron in false colour.
		
		*NB:* It is present when discretisation is available for the selected mesh generation method.  
	
	
	For physical cells, most entries have as many values as there are compartments.
	Spatial attributes (position coordinates, distance, and such) are given in microns following the NeuroML convention.
	
	For thinking of the neuron as a "tree" data structure, refer to:
	
	- the *Book of Genesis*, chapter 5, section 5: http://www.genesis-sim.org/GENESIS/iBoG/iBoGpdf/chapt5.pdf 
	- the TREES toolbox for analysing the structure of neurons: https://doi.org/10.1371/journal.pcbi.1000877
	- and the related literature.
	
	'''
	
	# TODO kwargs like selected cell names!
	# LATER allow exporting each cell type as an json and obj file ... or not?
	args = ["eden", "nml", nml_file, "explain", "cell", json.dumps(kwargs) ]
	# TODO use a temp file!
	out_file = 'explain_cell.json'
	_invoke_eden(args=args, verbose=verbose, threads=threads,
		extra_cmdline_args=extra_cmdline_args, executable_path=executable_path,  full_cmdline=full_cmdline)
	# load json file!
	with open(out_file,'r') as f:
		out_dict = json.load(f)
	
	# convert most lists to numpy nd_arrays
	for _, cell_info in out_dict.items():
		for key, val in cell_info.items(): # LATER for selected keys...
			if key in cell_info and key not in ['segment_groups']:
				cell_info[key] = np.array(val)
	
	return out_dict

def GetLemsLocatorsForCell(cell_info, compartment_ids=None):
	'''
	Generate a list of (extended) LEMS locators to access a property all over a cell.

	Parameters
	-----
	cell_info: dict

		The discretisation information about the cell, as returned from `eden_simulator.experimental.explain_cell`.
	
	compartment_ids: iterable
	
		A specific list of compartments to access (by default all over the cell)
	
	Returns
	-------
	l : list(str)

		The list of locators, to be used in extended LEMS paths capturing the cell property.
	
	See Also
	---
	:py:mod:`eden_simulator.experimental.explain_cell`: For obtaining ``cell_info``.
	`LEMS paths for cell locations <extension_paths.rst#lems-paths-for-cell-locations>`__ : For LEMS segment locators including `fractionAlong`.
	:doc:`intro_spatial`: To lean more about spatially detailed cells.
	
	'''
	comp_mid_seg, comp_mid_fra = (cell_info[x] for x in ('comp_midpoint_segment', 'comp_midpoint_fractionAlong'))
	if compartment_ids is None: compartment_ids = range(len(comp_mid_seg))
	return [str(comp_mid_seg[i])+('%.9f'%comp_mid_fra[i])[1:] for i in compartment_ids]
