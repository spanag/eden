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

def runEden( example_lems_file, threads=None, extra_cmdline_args=None, executable_path=None, verbose = False, full_cmdline=None):
	'''
	Run LEMS/NeuroML2 file with EDEN
	'''
	cwd = os.getcwd()
	#print(cwd)
	
	my_env = os.environ.copy()
	# my_env["LD_LIBRARY_PATH"] = cwd # LATER add this to the actual location of the executable, in case it is not delocated enough...
	if "OMP_SCHEDULE" not in my_env:
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
		# more places to check: https://github.com/PyTables/PyTables/blob/v3.8.0/tables/utils.py#L402 
		
		if threads > 8:
			# some setups, some times, have horrible slowdown when OMP_NUM_THREADS == logical cores, go figure
			# excluding one logical core shouldn't impact performance much due to the memory wall
			threads = threads - 1
			
	my_env["OMP_NUM_THREADS"] = str(threads)
	
	args = ["eden", "nml", example_lems_file, "gcc"]
	
	# If Eden executable is bundled in the package, use that by default 
	eden_bundled_exe_filename = embeden.get_exe_path()
	
	# But in all cases, executable_path takes top priority
	if executable_path:
		args[0] = executable_path
		
	elif eden_bundled_exe_filename:
		if verbose:
			print("Using bundled executable: "+eden_bundled_exe+" on "+eden_bundled_exe_filename)
		args[0] = eden_bundled_exe_filename
	
	if extra_cmdline_args:
		args += extra_cmdline_args
	if full_cmdline:
		args = full_cmdline
	
	if verbose:
		print(args)
	
	tic = time.time()
	t_run = datetime.now()
	
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
	
	results_Eden = reload_saved_data(example_lems_file, t_run=t_run)
	
	return OrderedDict(sorted(results_Eden.items()))


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

    if reload_events:
        event_output_files = sim.findall(ns_prefix + "EventOutputFile")
        for i, of in enumerate(event_output_files):
            name = of.attrib["fileName"]
            
            file_name = os.path.join(base_dir, name)
            if not os.path.isfile(file_name):  # If not relative to the LEMS file...
                file_name = os.path.join(base_lems_file_path, name)
            # if not os.path.isfile(file_name): # If not relative to the LEMS file...
            #    file_name = os.path.join(os.getcwd(),name)
            # ... try relative to cwd.
            # if not os.path.isfile(file_name): # If not relative to the LEMS file...
            #    file_name = os.path.join(os.getcwd(),'NeuroML2','results',name)
            # ... try relative to cwd in NeuroML2/results subdir.
            
            if not os.path.isfile(file_name):  # If not relative to the base dir...
                raise OSError(
                    ("Could not find simulation output event " "file %s" % file_name)
                )
            format = of.attrib["format"]
            logger.info("Loading saved events from %s (format: %s)" % (file_name, format))
            selections = {}
            for col in of.findall(ns_prefix + "EventSelection"):
                id = int(col.attrib["id"])
                select = col.attrib["select"]
                events[select] = []
                selections[id] = select

            with open(file_name) as f:
                for line in f:
                    values = line.split()
                    # TODO improve procssing speed
                    if format == "TIME_ID":
                        t = float(values[0])
                        id = int(values[1])
                    elif format == "ID_TIME":
                        id = int(values[0])
                        t = float(values[1])
                #logger.debug("Found a event in cell %s (%s) at t = %s" % (id, selections[id], t))
                events[selections[id]].append(t)

            if remove_dat_files_after_load:
                logger.warning( "Removing file %s after having loading its data!" % file_name)
                os.remove(file_name)

    output_files = sim.findall(ns_prefix + "OutputFile")
    n_output_files = len(output_files)

    for i, of in enumerate(output_files):
        traces["t"] = []
        name = of.attrib["fileName"]
        
        file_name = os.path.join(base_dir, name)
        if not os.path.isfile(file_name):  # If not relative to the LEMS file...
            file_name = os.path.join(base_lems_file_path, name)
        if not os.path.isfile(file_name):  # If not relative to the LEMS file...
            file_name = os.path.join(os.getcwd(), name) # ... try relative to cwd.
        if not os.path.isfile(file_name):  # If not relative to the LEMS file...
            file_name = os.path.join(os.getcwd(), "NeuroML2", "results", name) # ... try relative to cwd in NeuroML2/results subdir.
        
        if not os.path.isfile(file_name):  # If not relative to the LEMS file...
            raise OSError(("Could not find simulation output " "file %s" % file_name))
        
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
