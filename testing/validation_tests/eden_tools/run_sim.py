# Core Libraries
from collections import OrderedDict

# System libraries
import os
import os.path

import time

# PyNeuroML
from pyneuroml import pynml


def runNeuron( example_lems_file, verbose=False, reload_events=False):
	'''
	Convert LEMS/NeuroML2 file to NEURON with jNeuroML & run
	'''
	# Beware: jNeuroML may output files along with .xml file, instead of cwd like EDEN does.
	# reload_saved_data could get confused with the possible locations, pray it doesn't !
	out_dir,rel_filename = os.path.split(example_lems_file)
	
	tic = time.time()
	results = pynml.run_lems_with_jneuroml_neuron(rel_filename, nogui=True, load_saved_data=True, reload_events=reload_events, exec_in_dir=out_dir)
	if reload_events: traje, event = results # decompose tuple
	else: traje, event = (results, None)
	if not (traje or event):
		raise RuntimeError('Could not run simulation')
	toc = time.time()
	if verbose:
		print( "Ran jNeuroML_NEURON in %.2f seconds" % (toc - tic) )
	
	if reload_events:
		return OrderedDict(sorted(traje.items())), OrderedDict(sorted(event.items()))
	else:
		return OrderedDict(sorted(traje.items()))

def runJLems( example_lems_file, verbose=False, reload_events=False):
	'''
	Run LEMS/NeuroML2 file with jLEMS
	'''
	# Beware: jNeuroML may output files along with .xml file, instead of cwd like EDEN does.
	# reload_saved_data could get confused with the possible locations, pray it doesn't !
	out_dir,rel_filename = os.path.split(example_lems_file)
	
	tic = time.time()
	results = pynml.run_lems_with_jneuroml(rel_filename, nogui=True, load_saved_data=True, reload_events=reload_events, exec_in_dir=out_dir)
	if reload_events: traje, event = results # decompose tuple
	else: traje, event = (results, None)
	if not (traje or event):
		raise RuntimeError('Could not run simulation')
	toc = time.time()
	if verbose:
		print( "Ran jNeuroML_LEMS in %.2f seconds" % (toc - tic) )
	if reload_events:
		return OrderedDict(sorted(traje.items())), OrderedDict(sorted(event.items()))
	else:
		return OrderedDict(sorted(traje.items()))
