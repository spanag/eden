# Core Libraries
from collections import OrderedDict

# System libraries
import subprocess
import os
import os.path
	
import time
from datetime import datetime

# NumPy
import numpy as np

# Some logging built in Python
import logging
logger = logging.getLogger(__name__)

# Own files
from . import embeden

def _autodetect_threads():
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
	# more places to check: https://github.com/PyTables/PyTables/blob/v3.8.0/tables/utils.py#L402 
	
	# binder and other environments may allocate 1 thread or less out of the many in the visible cpu: https://discourse.jupyter.org/t/mybinder-and-multiprocessing/3238
	# Trying to multi-thread in each step thus becomes near unusable.
	# TODO auto-detect the situation and handle in Eden proper.
	if 'BINDER_LAUNCH_HOST' in os.environ:
		threads = 1
	if 'READTHEDOCS' in os.environ:
		threads = 1 # just in case
	if 'DEEPNOTE_CPU_COUNT' in os.environ:
		deepcpus = os.environ['DEEPNOTE_CPU_COUNT']
		if deepcpus.isdigit(): threads = int(deepcpus)
	
	if threads > 8:
		# some setups, some times, have horrible slowdown when OMP_NUM_THREADS == logical cores, go figure
		# excluding one logical core shouldn't impact performance much due to the memory wall
		threads = threads - 1
	return threads

def _invoke_eden(args, threads=None, extra_cmdline_args=None, executable_path=None, verbose = False, full_cmdline=None):
	cwd = os.getcwd()
	# print(cwd)
	args = list(args) # clone
	my_env = os.environ.copy()
	# my_env["LD_LIBRARY_PATH"] = cwd # LATER add this to the actual location of the executable, in case it is not delocated enough...
	if "OMP_SCHEDULE" not in my_env:
		my_env["OMP_SCHEDULE"] = "static"
	
	if threads is None:
		threads = _autodetect_threads()
			
	my_env["OMP_NUM_THREADS"] = str(threads)
	
	
	# If Eden executable is bundled in the package, use that by default 
	eden_bundled_exe_filename = embeden.get_exe_path()
	
	# But in all cases, executable_path takes top priority
	if executable_path:
		args[0] = executable_path
		
	elif eden_bundled_exe_filename:
		if verbose:
			print("Using bundled executable: on "+eden_bundled_exe_filename)
		args[0] = eden_bundled_exe_filename
	
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
	

def runEden( example_lems_file, *,
	reload_events=False, verbose = False, threads = None,
	extra_cmdline_args=None, full_cmdline=None, executable_path=None
):
	'''
	Run a NeuroML2/LEMS file with EDEN.
	
	
	Parameters
	----------
	
	example_lems_file : str
		The filename where the `<Simulation>` to run is located.
	
	reload_events : bool, optional
		Whether to return a:
			- `dict` with the recorded time series, mapping each recorded `<LemsQuantityPath>` to its sampled values over time. All samples are taken at the same points in time, represented by the additional time vector `t` in seconds.
			
			- or a pair of `dict` s, the first as before and the second mapping each recorded event to its times of occurrence.
	
	threads : int, optional
		The number of threads to run the simulation with, or `None` for automatic selection. 
	
	verbose : bool, optional
		Add some more console output when running the simulator.
	
	extra_cmdline_args : list[str], optional
		Additional command line arguments to pass to EDEN.  
		Refer to the EDEN user's manual, "Command line API" for more details.
		
	full_cmdline : list[str], optional
		Specify the exact argv to be passed to EDEN.  
		Refer to the EDEN user's manual, "Command line API" for more details.
	
	executable_path : str, optional
		Select a specific executable of EDEN to run the simulation with. Useful for custom installations and uses of EDEN.
	
	
	Returns
	-------

	trajectories : dict[str, numpy.ndarray[float[n]]]
		
		- *when reload_events == False*
		
	(trajectories, events) :  tuple(dict, dict)
		
		- *when reload_events == True*; see `reload_events` parameter.
	
	
	Notes
	-----
	
	All values are returned in `fundamental` units or derived thereof.  
	These are not always SI units!
	
	- e.g. molarity is in `mol/m3`, not `mol/L` !
		
	Refer to https://en.wikipedia.org/wiki/SI_derived_unit for a list of such units.
	'''
	
	args = ["eden", "nml", example_lems_file, "gcc"]
	
	t_run = datetime.now()
	_invoke_eden(args, threads, extra_cmdline_args, executable_path, verbose, full_cmdline)
	
	results = reload_saved_data(example_lems_file, reload_events=reload_events, t_run=t_run)
	if reload_events: traje, event = results # decompose tuple
	else: traje, event = (results, None)
	if reload_events:
		return OrderedDict(sorted(traje.items())), OrderedDict(sorted(event.items()))
	else:
		return OrderedDict(sorted(traje.items()))


# Adapted from pynml.reload_saved_data: https://github.com/NeuroML/pyNeuroML/blob/v0.7.5/pyneuroml/pynml.py
def reload_saved_data(
	lems_file_name,
	base_dir=".",
	t_run=datetime(1900, 1, 1),
	reload_events=False,
	remove_dat_files_after_load=False,
):
	"""Reload data saved from previous EDEN simulation run.
	
	:param lems_file_name: name of LEMS file that was used to run the simulation.
	:type lems_file_name: str
	:param base_dir: directory to run in
	:type base_dir: str
	:param t_run: time of run
	:type t_run: datetime
	:param reload_event: toggle whether events should be loaded
	:type reload_event: bool
	:param remove_dat_files_after_load: toggle if data files should be deleted after they've been loaded
	:type remove_dat_files_after_load: bool
	"""
	
	if not os.path.isfile(lems_file_name):
		real_lems_file = os.path.realpath(os.path.join(base_dir, lems_file_name))
	else:
		real_lems_file = os.path.realpath(lems_file_name)

	logger.debug(
		"Reloading data specified in LEMS file: %s (%s), base_dir: %s, cwd: %s"
		% (lems_file_name, real_lems_file, base_dir, os.getcwd())
	)
	# Could use pylems to parse all this...
	from collections import OrderedDict
	traces = OrderedDict()
	events = OrderedDict()

	base_lems_file_path = os.path.dirname(os.path.realpath(lems_file_name))
	from lxml import etree
	tree = etree.parse(real_lems_file)

	sim = tree.getroot().find("Simulation")
	ns_prefix = ""

	possible_prefixes = ["{http://www.neuroml.org/lems/0.7.2}", "{http://www.neuroml.org/schema/neuroml2}"]
	if sim is None:
		for pre in possible_prefixes:
			for comp in tree.getroot().findall(pre + "Component"):
				if comp.attrib["type"] == "Simulation":
					ns_prefix = pre
					sim = comp
	
	def FindAllFromMulti(root_tag, tagnames):
		return [ x for tagname in tagnames for x in list(root_tag.findall(ns_prefix + tagname))  ]
	
	def GetFileNameLocationFromOutputTag(of):
		try:
			name = of.attrib["fileName"]
		except KeyError:
			name = of.attrib["href"]
			if name.startswith('file://'):
				name = name[len('file://'):] # just chop it off
			elif name.find('://') >= 0:
				return None # ignore, for it could be a funky stream. Just don't reload it.
			else: name = base_lems_file_path+'/'+name# For Eden's "href" format, no URI scheme means: file relative to sim file !
				# raise ValueError('Cannot reload path '+name) # unknown url scheme
		file_name = os.path.join(base_dir, name)
		return file_name
	
	if reload_events:
		event_output_files = FindAllFromMulti(sim, ["EventOutputFile","EdenEventOutputFile"])
		for i, of in enumerate(event_output_files):

			file_name = GetFileNameLocationFromOutputTag(of)
			if not os.path.isfile(file_name):  # If not relative to the base dir...
				continue # ignore, for it could be a funky stream. Just don't reload it.
				# TODO discern whether it is a magical file specifically, and complain if it's totally missing, or a folder...
				# raise OSError(("Could not find simulation output event " "file %s" % file_name))
			
			format = of.attrib["format"]
			logger.info("Loading saved events from %s (format: %s)" % (file_name, format))
			selections = {}
			for col in of.findall(ns_prefix + "EventSelection"):
				id = int(col.attrib["id"])
				select = col.attrib["select"]
				events[select] = []
				selections[id] = select
			
			# FIXME eden's formats as well !!
			with open(file_name) as f:
				for line in f:
					values = line.split()
					# TODO improve processing speed
					if format == "TIME_ID":
						t = float(values[0])
						id = int(values[1])
					elif format == "ID_TIME":
						id = int(values[0])
						t = float(values[1])
					else: raise ValueError(format)
					#logger.debug("Found a event in cell %s (%s) at t = %s" % (id, selections[id], t))
					events[selections[id]].append(t)

			if remove_dat_files_after_load:
				logger.warning( "Removing file %s after having loading its data!" % file_name)
				os.remove(file_name)

	output_files = FindAllFromMulti(sim, ["OutputFile","EdenOutputFile"])
	n_output_files = len(output_files)

	for i, of in enumerate(output_files):
		
		file_name = GetFileNameLocationFromOutputTag(of)
		if not os.path.isfile(file_name):  # If not relative to the LEMS file...
			continue # ignore, for it could be a funky stream. Just don't reload it.
			# raise OSError(("Could not find simulation output " "file %s" % file_name))
		
		t_file_mod = datetime.fromtimestamp(os.path.getmtime(file_name))
		if t_file_mod < t_run:
			raise Exception(
				"Expected output file %s has not been modified since "
				"%s but the simulation was run later at %s."
				% (file_name, t_file_mod, t_run)
			)

		logger.debug("Loading saved data from %s" % (file_name))

		cols = []
		cols.append("t")
		for col in of.findall(ns_prefix + "OutputColumn"):
			quantity = col.attrib["quantity"]
			traces[quantity] = []
			cols.append(quantity)

		with open(file_name) as f:
			data = np.fromstring(f.read(),sep=' ').reshape((-1, len(cols)), order='C')
			for i,col in enumerate(cols):
				traces[col] = data[:,i]
			
		if remove_dat_files_after_load:
			logger.warning("Removing file %s after having loading its data!" % file_name)
			os.remove(file_name)
	
	if reload_events:
		return traces, events
	else:
		return traces
