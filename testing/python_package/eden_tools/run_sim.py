# Core Libraries
from collections import OrderedDict

# System libraries
import sys
import subprocess
import os
import os.path

import time

# PyNeuroML
from pyneuroml import pynml


def runNeuron( example_lems_file, verbose=False):
	'''
	Convert LEMS/NeuroML2 file to NEURON with jNeuroML & run
	'''
	# Beware: jNeuroML may output files along with .xml file, instead of cwd like EDEN does.
	# reload_saved_data could get confused with the possible locations, pray it doesn't !
	out_dir,rel_filename = os.path.split(example_lems_file)
	
	tic = time.time()
	results_Neuron = pynml.run_lems_with_jneuroml_neuron(rel_filename, nogui=True, load_saved_data=True, exec_in_dir=out_dir)
	if not results_Neuron:
		raise RuntimeError('Could not run simulation')
	toc = time.time()
	if verbose:
		print( "Ran jNeuroML_NEURON in %.2f seconds" % (toc - tic) )
	
	return OrderedDict(sorted(results_Neuron.items()))

def runJLems( example_lems_file, verbose=False):
	'''
	Run LEMS/NeuroML2 file with jLEMS
	'''
	# Beware: jNeuroML may output files along with .xml file, instead of cwd like EDEN does.
	# reload_saved_data could get confused with the possible locations, pray it doesn't !
	out_dir,rel_filename = os.path.split(example_lems_file)
	
	tic = time.time()
	results = pynml.run_lems_with_jneuroml(rel_filename, nogui=True, load_saved_data=True, exec_in_dir=out_dir)
	if not results:
		raise RuntimeError('Could not run simulation')
	toc = time.time()
	if verbose:
		print( "Ran jNeuroML_LEMS in %.2f seconds" % (toc - tic) )
	
	return OrderedDict(sorted(results.items()))

def runEden( example_lems_file, threads=None, extra_cmdline_args=None, executable_path=None, verbose = False, full_cmdline=None):
	'''
	Run LEMS/NeuroML2 file with EDEN
	'''
	cwd = os.getcwd()
	#print(cwd)
	
	my_env = os.environ.copy()
	# my_env["LD_LIBRARY_PATH"] = cwd
	my_env["OMP_SCHEDULE"] = "static"
	
	if threads is None:
		
		try: 
			import multiprocessing
		except ImportError:
			multiprocessing = None
		
		if hasattr(os, 'sched_getaffinity'):
			threads = len(os.sched_getaffinity(0)) # XXX still just a guess, also doesn't take into account the amount of physical cores, etc.
		elif multiprocessing is not None and hasattr(multiprocessing, 'cpu_count'):
			threads = multiprocessing.cpu_count()
		else:
			# multiprocessing.cpu_count is missing, this is strange
			# Select one thread, let user override
			threads = 1
		
		if threads > 8:
			# some setups, some times, have horrible slowdown when OMP_NUM_THREADS == logical cores, go figure
			# excluding one logical core shouldn't impact performance much due to the memory wall
			threads = threads - 1
			
	my_env["OMP_NUM_THREADS"] = str(threads)
	
	args = ["eden", "nml", example_lems_file, "gcc"]
	
	# If Eden executable is bundled in the package, use that by default 
	import pkg_resources
	import platform
	exe_extension = ".exe" if platform.system() == 'Windows' else ""
	eden_bundled_exe = "../bin/eden"+exe_extension
	
	if pkg_resources.resource_exists(__name__, eden_bundled_exe):
		eden_bundled_exe_filename = pkg_resources.resource_filename(__name__, eden_bundled_exe)
		if verbose:
			print("Using bundled executable: "+eden_bundled_exe+" on "+eden_bundled_exe_filename)
		args[0] = eden_bundled_exe_filename
	
	if executable_path:
		args[0] = executable_path
	
	if extra_cmdline_args:
		args += extra_cmdline_args
	if full_cmdline:
		args = full_cmdline
	
	if verbose:
		print(args)
	
	tic = time.time()
	
	sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env = my_env)
	out, err = map( lambda out: out.decode('utf-8'), sp.communicate() )
	if sp.returncode != 0 :
		print("returncode of subprocess:")
		print(sp.returncode)
		if out:
			print("standard output of subprocess:")
			print(out)
		if err:
			print("standard error of subprocess:")
			print(err)
		raise RuntimeError('Could not run simulation')

	toc = time.time()
	if verbose:
		print( "Ran EDEN in %.2f seconds" % (toc - tic) )
	
	results_Eden = pynml.reload_saved_data(example_lems_file, simulator = 'jNeuroML_NEURON');
	
	return OrderedDict(sorted(results_Eden.items()))
