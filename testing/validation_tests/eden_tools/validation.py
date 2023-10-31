import numpy as np

from eden_simulator import runEden
from .run_sim import runNeuron

# TODO add docstrings
def GetOverUnder_Box( results_Truth, time_Truth, dt_box, dv_box ):
	import scipy
	from scipy import ndimage
	
	dt = time_Truth[2] - time_Truth[1]
	# dt_box = 0.200 * 0.001
	# dv_box = 5 * 0.001
	tol_samples = dt_box / dt
	# print(tol_samples)
	strel_width = int(tol_samples) + 1
	
	dmax = scipy.ndimage.filters.maximum_filter1d( results_Truth, strel_width, mode='nearest') + dv_box
	dmin = scipy.ndimage.filters.minimum_filter1d( results_Truth, strel_width, mode='nearest') - dv_box
	
	return(dmax, dmin)
	
def VerifyTimeSeries_Box( results_Truth, results_Test, time_Truth, dt_box, dv_box ):
	'''Returns True if test passed; string describing first fault otherwise'''    
	results_Truth
	# XXX this won't work with irregular timesteps!
	# plus it allows for an extra time-error of dt due to rounding
	# TODO add a hand-written, sliding-window loop that also interpolates on the exact points the window ends
	# also check https://stackoverflow.com/a/43288787 strided arrays
	# and https://web.archive.org/web/20130615115546/http://www.richardhartersworld.com/cri/2001/slidingmin.html O(n) algorithm
	
	dmax, dmin = GetOverUnder_Box( results_Truth, time_Truth, dt_box, dv_box )
	
	over_points  = np.flatnonzero( results_Test > dmax )
	under_points = np.flatnonzero( results_Test < dmin )
	nan_points   = np.flatnonzero( np.isnan(results_Truth) | np.isnan(results_Test) )
	# print (over_points)
	# print (under_points)
	# print (nan_points)
	over_earliest = len(dmax) + 5 # invalid value, later than possible
	if over_points.size > 0:
		over_earliest = over_points[0]

	under_earliest = len(dmin) + 5 # invalid value, later than possible
	if under_points.size > 0:
		under_earliest = under_points[0]
	
	nan_earliest = len(dmin) + 5 # invalid value, later than possible
	if nan_points.size > 0: nan_earliest = nan_points[0]
	# print(over_earliest)
	if over_points.size > 0 or under_points.size > 0 or nan_points.size > 0:
		# one of them should be in the sample range
		if over_earliest < under_earliest and over_earliest < nan_earliest:
			fail_desc = "Over limit"
			fail_sample = over_earliest
			diff = (results_Test[over_earliest] - results_Truth[over_earliest])
		elif under_earliest < over_earliest and under_earliest < nan_earliest:
			fail_desc = "Under limit"
			fail_sample = under_earliest
			diff = (results_Test[under_earliest] - results_Truth[under_earliest])
		elif nan_earliest < over_earliest and nan_earliest < under_earliest:
			fail_desc = "NaN"
			fail_sample = nan_earliest
			diff = (results_Test[nan_earliest] - results_Truth[nan_earliest])
		else:
			raise AssertionError( 'exceeding both max and min limits')
		return {'desc': fail_desc, 'sample':fail_sample, 'time': time_Truth[fail_sample], 'diff':diff}
	
	else:
		return True
		
def VerifySpikeTrain_SameCount_Delay( results_Truth, results_Test, dt_box ):
	'''Returns True if test passed; string describing first fault otherwise'''
		
	if len(results_Truth) != len(results_Test):
		return {'desc': 'Wrong count', 'time': 0, 'diff': (len(results_Truth) - len(results_Test))  } # TODO more accurate time...
	
	results_Truth = np.array(results_Truth); results_Test = np.array(results_Test)
	dif = results_Test - results_Truth
	check = abs(dif) > dt_box
	if any(check):
		first_point = np.nonzero(check)[0][0]
		d = dif[first_point]
		# print(first_point,d)
		return {'desc': ("Before" if d < 0 else 'After')+" limit", 'time': results_Test[first_point], 'diff':d}
	else:
		return True

def ExplainVerification( Truth, Test, criteria, key ):
	
	from matplotlib import pyplot as plt
	
	results_Truth = Truth[key]
	results_Test = Test[key]
	
	time_Truth = Truth['t']
	
	if criteria == 'exact' :
		criterion = { 'type': 'box', 'dt': 0.0, 'dv': 0.0 }
	elif key in criteria:
		criterion = criteria[key]
	else:
		raise ValueError('key not found in criteria')

	typ = criterion['type']
	
	
	if typ == 'box' or typ == 'exact':
		if typ == 'box': dt, dv = criterion['dt'], criterion['dv']
		if typ == 'exact': dt, dv = 0, 0
		
		dmax, dmin = GetOverUnder_Box( results_Truth, time_Truth, dt,dv )
		# TODO the Truth in general, there might be another source of Truth (ikr)
		plt.plot( time_Truth, results_Truth, label="NEURON")
		plt.plot( time_Truth, results_Test, label="EDEN")
		plt.plot( time_Truth, dmax, linestyle='--', label='Max')
		plt.plot( time_Truth, dmin, linestyle='--', label='Min')
		res = VerifyTimeSeries_Box( results_Truth, results_Test, time_Truth, dt, dv )
		if res is not True:
			plt.axvline( res['time'], color='red', linewidth='1')
			plt.axhline( results_Test[res['sample']], color='red', linewidth='1')
	else:
		raise ValueError('Unknown criterion type "'+ typ + '"')
			
	
	plt.legend()

def VerifySimResults(
	results_Neuron, results_Eden, criteria,
	events_Neuron={}, events_Eden={},
	skip_missing_criteria=False,
	explain_verification_on_fail=False, 
):
	trajes = results_Eden.keys()
	time_Neuron = results_Neuron['t']
	time_Eden   = results_Eden  ['t']
	
	events = events_Eden.keys()
	
	ok = True
	
	validation_fail_trajes = []
	validation_fail_events = []
	validation_missing = []
	
	for key in trajes:
		if key == 't':
			continue
		res_neuron = results_Neuron[key]
		res_interp = np.interp(time_Neuron, time_Eden, results_Eden[key])
		
		# diff = np.subtract(res_neuron, res_interp)
		# if np.mean(np.abs(diff)) > 0.00001:
		# 	plt.plot(time_Neuron, diff, label="Diff: "+key)
		res = None
		criterion = None
		if criteria == 'exact' :
			criterion = { 'type': 'box', 'dt': 0.0, 'dv': 0.0 }
		elif key in criteria:
			criterion = criteria[key]
		else:
			if skip_missing_criteria:
				criterion = None
				pass
			else:
				# print("Validation test missing for "+ key)
				validation_missing.append(key)
		
		if criterion:
			if isinstance(criterion, str): typ = criterion
			else: typ = criterion['type']
			
			if typ == 'box':
				res = VerifyTimeSeries_Box( res_neuron, res_interp, time_Neuron, criterion['dt'], criterion['dv'])
			elif typ == 'exact':
				res = VerifyTimeSeries_Box( res_neuron, res_interp, time_Neuron, 0.0, 0.0)
			else:
				raise ValueError('Unknown criterion type "'+ typ + '"')
			VerifySpikeTrain_SameCount_Delay
			if res is None:
				# indeterminate??
				raise
			
			if res is not True:
				# it failed
				failed = dict(res)
				failed['key'] = key
				validation_fail_trajes.append(failed)
	
	for key in events:
		no_tolerance = 0e-9 # 1 nsec to avoid printf roundoff? TODO investigate
		res_neuron = events_Neuron[key]
		res__eden = events_Eden[key]
		# print(key, events)
		exact = False
		criterion = None
		if criteria == 'exact' :
			criterion = 'exact'
		elif key in criteria:
			criterion = criteria[key]
		else:
			if skip_missing_criteria:
				criterion = None
				pass
			else:
				# print("Validation test missing for "+ key)
				validation_missing.append(key)
		# print(criterion, criteria, key in criteria)
		if criterion:
			if isinstance(criterion, str): typ = criterion
			else: typ = criterion['type']
			
			if typ == 'exact':
				res = VerifySpikeTrain_SameCount_Delay(res_neuron, res__eden, no_tolerance)
			elif typ == 'samecount_box':
				res = VerifySpikeTrain_SameCount_Delay(res_neuron, res__eden, criterion['dt'])
			else:
				raise ValueError('Unknown criterion type "'+ typ + '"')
			
			if res is None:
				# indeterminate??
				raise
			
			if res is not True:
				# it failed
				failed = dict(res)
				failed['key'] = key
				validation_fail_events.append(failed)
	
	if validation_missing:
		ok = False
		print("Missing time series validation criteria:")
		
		for key in validation_missing:
			print(f"\t{key}\n")
	
	validation_fail_trajes.sort(key=lambda x: x['time'])
	validation_fail_events.sort(key=lambda x: x['time'])
	
	if validation_fail_trajes:
		ok = False

		print("Failed time series:")

		for res in validation_fail_trajes:

			resstr = res['key'] + ' : ' + res['desc']  + " at "+ ("%.3f msec" % (res['time'] * 1000) ) + (" (sample %d, diff %f)" % (res['sample'], res['diff']) )
			print("\t" + resstr)

		if explain_verification_on_fail:
			failed = validation_fail_trajes[0]
			ExplainVerification(results_Neuron, results_Eden, criteria, failed['key'])
			from matplotlib import pyplot as plt
			plt.show()
			
	if validation_fail_events:
		ok = False
		print("Failed spike trains:")

		for res in validation_fail_events:

			resstr = res['key'] + ' : ' + res['desc']  + " at "+ ("%.3f msec" % (res['time'] * 1000) ) + (" (diff %f)" % (res['diff']) )
			print("\t" + resstr)

		if explain_verification_on_fail:
			failed = validation_fail_events[0]
			# ExplainVerification(results_Neuron, results_Eden, criteria, failed['key'])
			from matplotlib import pyplot as plt; plt.show()
		
	
	return ok

def RunTest(
	test,
	skip_missing_criteria = False,
	explain_verification_on_fail=False,
	extra_cmdline_args=None,
	verbose = False,
	executable_path=None,
):
	ok = False # return value
	# since this is testing code,
	# for debugging the test suite without rebuilding, slip eden/testing/sandbox into the PATH
	import sys
	if True and any(('testing/sandbox' in x) for x in sys.path):	
		from importlib import reload;
		from wheel_hollow.eden_simulator import run_sim; reload(run_sim); from wheel_hollow.eden_simulator.run_sim import runEden
		import eden_tools.run_sim as run_sim; reload(run_sim); from eden_tools.run_sim import runNeuron
	else:
		global runEden, runNeuron
	
	criteria = test.get('validation_criteria', None)
	
	skip_missing = skip_missing_criteria
	if skip_missing_criteria is not True:
		skip_missing = not not test.get('skip_missing',False)
	
	typ = test['type']
	if typ == 'smoke_test':
		test_filename = test['sim_file']
		print('Checking if %s runs...' % test_filename)
		try:
			res, eve = runEden( test_filename, reload_events=True, executable_path=executable_path )
			ok = True
			print('PASS %s smoke test' % test_filename)
		except RuntimeError:
			print('FAIL %s smoke test' % test_filename)
		
	elif(
		typ == 'validate_vs_neuron'
		or typ == 'eden_vs_eden'
		or typ == 'half_vs_half'
		or typ == 'function_vs_function'
	):
		# TODO add jLEMS
		test_filename = test['sim_file']
		
		testcase_type = 'test'
		if typ == 'validate_vs_neuron':
			testcase_type = 'vs. NEURON'
			print('Validating %s %s...' % ( test_filename, testcase_type ) )
			results_Neuron, events_Neuron = runNeuron( test_filename, verbose=verbose, reload_events=True )
			
			eden_extra_kwargs = test.get('eden_extra_kwargs', {})
			results_Eden  , events_Eden   = runEden  ( test_filename, verbose=verbose, extra_cmdline_args=extra_cmdline_args, reload_events=True, executable_path=executable_path  , **eden_extra_kwargs )
			
		elif typ == 'eden_vs_eden':
			testcase_type = 'vs. EDEN'
			
			default_kwargs = { 'verbose': verbose, 'extra_cmdline_args': extra_cmdline_args }
			
			truth_kwargs = test.get('truth_kwargs',default_kwargs)
			test_kwargs  = test.get('test_kwargs', default_kwargs)
			
			print( 'Validating %s, one run of EDEN (%s) with another(%s)...' % ( test_filename, str(truth_kwargs), str(test_kwargs) ) )
			results_Neuron, events_Neuron = runEden( test_filename, reload_events=True, executable_path=executable_path , **truth_kwargs )
			results_Eden  , events_Eden   = runEden( test_filename, reload_events=True, executable_path=executable_path , **test_kwargs )
			
		elif typ == 'half_vs_half':
			testcase_type = 'vs. Half'
			test_filename = test['sim_file']
			compare = test['validation_criteria']
			ignore = test.get('ignore',[])
			for pair in compare:
				if not (
					isinstance(pair, list)
					and len(pair) == 3
					and isinstance(pair[0], str)
					and isinstance(pair[1], str)
					and ( isinstance(pair[2], dict) or pair[2] == 'exact')
				):
					raise ValueError('elements of half_vs_half "validation_criteria" must be triplets of string series truth, string series test, dict comparison criterion')
			
			default_kwargs = { 'verbose': verbose, 'extra_cmdline_args': extra_cmdline_args }
			test_kwargs = default_kwargs
			
			print( 'Validating %s, trajectories in one run of EDEN (%s)' % ( test_filename, str(test_kwargs) ) )
			
			results, events  = runEden( test_filename, reload_events=True, executable_path=executable_path , **test_kwargs )
			# print(events)
			results_test = {'t':results['t']}
			results_true = {'t':results['t']}
			events__true = {}
			events__test = {}
			results_tezt = {} # just to get the set of trajectories in the 'test' set :D
			criteria = {}
			remainder = {}
			
			for pair in compare:
				ltrue = pair[0]
				ltest = pair[1]
				if ltrue in results:
					results_true[ltrue] = results[ltrue]
					results_test[ltrue] = results[ltest]
					results_tezt[ltest] = results[ltest]
				elif ltrue in events:
					events__true[ltrue] = events[ltrue]
					events__test[ltrue] = events[ltest]
					results_tezt[ltest] = events[ltest]
				else:
					raise ValueError(f'log {ltrue} not found in trajes or events!')
				criteria[ltrue] = pair[2]
			
			# check for any remaining trajectories
			for key,val in list(results.items()) + list(events.items()):
				if key == "t": continue
				if key not in results_true and key not in results_tezt and key not in ignore:
					remainder[key] = val
			if not skip_missing_criteria and remainder:
				ok = False
				print("Some time series were not checked:")
				for key in remainder:
					print("\t" + key)
			
			default_kwargs = { 'verbose': verbose, 'extra_cmdline_args': extra_cmdline_args }
			
			results_Neuron = results_true; events_Neuron = events__true
			results_Eden   = results_test; events_Eden   = events__test
			
		elif typ == 'function_vs_function':
			testcase_type = 'vs. custom'
			
			default_kwargs = { 'verbose': verbose, 'extra_cmdline_args': extra_cmdline_args }
			
			results_Neuron = test['truth'](**default_kwargs)
			results_Eden   = test['test' ](**default_kwargs)
			#TODO use if needed...
		else:
			raise ValueError('unknown subtest type '+ typ+'??')
		
		# print(events_Neuron)
		res = VerifySimResults( results_Neuron, results_Eden, criteria,
			events_Neuron=events_Neuron, events_Eden=events_Eden,
			skip_missing_criteria=skip_missing, explain_verification_on_fail=explain_verification_on_fail)
		if res is True:
			ok = True
			print('PASS %s %s' % ( test_filename, testcase_type ) )
	else:
		raise ValueError('Unknown test type "'+ typ + '"')
	
	return ok

