import numpy as np

from .run_sim import runNeuron, runEden

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
	
	# XXX this won't work with irregular timesteps!
	# plus it allows for an extra time-error of dt due to rounding
	# TODO add a hand-written, sliding-window loop that also interpolates on the exact points the window ends
	# also check https://stackoverflow.com/a/43288787 strided arrays
	# and https://web.archive.org/web/20130615115546/http://www.richardhartersworld.com/cri/2001/slidingmin.html O(n) algorithm
	
	dmax, dmin = GetOverUnder_Box( results_Truth, time_Truth, dt_box, dv_box )
	
	over_points  = np.flatnonzero( results_Test > dmax )
	under_points = np.flatnonzero( results_Test < dmin )
	# print (over_points)
	# print (under_points)
	over_earliest = len(dmax) + 5 # invalid value, later than possible
	if over_points.size > 0:
		over_earliest = over_points[0]

	under_earliest = len(dmin) + 5 # invalid value, later than possible
	if under_points.size > 0:
		under_earliest = under_points[0]
	# print(over_earliest)
	# print(under_earliest)
	if over_points.size > 0 or under_points.size > 0:
		# one of them should be in the sample range
		if over_earliest < under_earliest:
			fail_desc = "Over limit"
			fail_sample = over_earliest
		elif under_earliest < over_earliest:
			fail_desc = "Under limit"
			fail_sample = under_earliest
		else:
			raise AssertionError( 'exceeding both max and min limits')
		return {'desc': fail_desc, 'sample':fail_sample, 'time': time_Truth[fail_sample]}
	
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

	if typ == 'box':
		dmax, dmin = GetOverUnder_Box( results_Truth, time_Truth, criterion['dt'], criterion['dv'] )
		# TODO the Truth in general, there might be another source of Truth (ikr)
		plt.plot( time_Truth, results_Truth, label="NEURON")
		plt.plot( time_Truth, results_Test, label="EDEN")
		plt.plot( time_Truth, dmax, linestyle='--', label='Max')
		plt.plot( time_Truth, dmin, linestyle='--', label='Min')
		res = VerifyTimeSeries_Box( results_Truth, results_Test, time_Truth, criterion['dt'], criterion['dv'] )
		if res is not True:
			plt.axvline( res['time'], color='red', linewidth='1')
			plt.axhline( results_Test[res['sample']], color='red', linewidth='1')
	else:
		raise ValueError('Unknown criterion type "'+ typ + '"')
			
	
	plt.legend()

def VerifySimResults(
	results_Neuron, results_Eden, criteria,
	skip_missing_criteria=False,
	explain_verification_on_fail=False, 
):
	plot_vars = results_Eden.keys()
	time_Neuron = results_Neuron['t']
	time_Eden   = results_Eden  ['t']
	
	
	ok = True
	
	validation_fails = []
	validation_missing = []
	
	for key in plot_vars:
		if key == 't':
			continue
		res_neuron = results_Neuron[key]
		res_interp = np.interp(time_Neuron, time_Eden, results_Eden[key])
		
#         diff = np.subtract(res_neuron, res_interp)
#         if np.mean(np.abs(diff)) > 0.00001:
#             plt.plot(time_Neuron, diff, label="Diff: "+key)
		res = None
		
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
			typ = criterion['type']
			
			if typ == 'box':
				res = VerifyTimeSeries_Box( res_neuron, res_interp, time_Neuron, criterion['dt'], criterion['dv'])
			else:
				raise ValueError('Unknown criterion type "'+ typ + '"')
			
			if res is None:
				# indeterminate??
				raise
			
			if res is not True:
				# it failed
				failed = dict(res)
				failed['key'] = key
				validation_fails.append(failed)
			
			
		
	if validation_missing:
		ok = False
		print("Missing time series validation criteria:")
		
		for key in validation_missing:
			print("\t" + key)
	
	validation_fails.sort(key=lambda x: x['time'])
	
	if validation_fails:
		ok = False

		print("Failed time series:")

		for res in validation_fails:

			resstr = res['key'] + ' : ' + res['desc']  + " at "+ ("%.3f msec" % (res['time'] * 1000) ) + (" (sample %d)" % res['sample'] )
			print("\t" + resstr)

		if explain_verification_on_fail:
			failed = validation_fails[0]
			ExplainVerification(results_Neuron, results_Eden, criteria, failed['key'])
			from matplotlib import pyplot as plt
			plt.show()
		
		
	return ok

def RunTest(
	test,
	skip_missing_criteria = False,
	explain_verification_on_fail=False,
	extra_cmdline_args=None,
	verbose = False,
):
	ok = False # return value
	
	typ = test['type']
	if typ == 'smoke_test':
		test_filename = test['sim_file']
		print('Checking if %s runs...' % test_filename)
		try:
			res = runEden( test_filename )
			ok = True
			print('PASS %s smoke test' % test_filename)
		except RuntimeError:
			print('FAIL %s smoke test' % test_filename)
		
	elif(
		typ == 'validate_vs_neuron'
		or typ == 'eden_vs_eden'
	):
		# TODO add jLEMS
		test_filename = test['sim_file']
		
		testcase_type = 'test'
		if typ == 'validate_vs_neuron':
			testcase_type = 'vs. NEURON'
			print('Validating %s %s...' % ( test_filename, testcase_type ) )
			results_Neuron = runNeuron( test_filename, verbose=verbose )
			results_Eden   = runEden  ( test_filename, verbose=verbose, extra_cmdline_args=extra_cmdline_args )
			
		elif typ == 'eden_vs_eden':
			testcase_type = 'vs. EDEN'
			
			default_kwargs = { 'verbose': verbose, 'extra_cmdline_args': extra_cmdline_args }
			
			truth_kwargs = test.get('truth_kwargs',default_kwargs)
			test_kwargs  = test.get('test_kwargs', default_kwargs)
			
			print( 'Validating %s, one run of EDEN (%s) with another(%s)...' % ( test_filename, str(truth_kwargs), str(test_kwargs) ) )
			results_Neuron = runEden( test_filename, **truth_kwargs )
			results_Eden   = runEden( test_filename, **test_kwargs  )
			
		else:
			raise ValueError('unknown subtest type '+ typ+'??')
		
		skip_missing = skip_missing_criteria
		if skip_missing_criteria is not True:
			skip_missing = not not test.get('skip_missing',False)
		criteria = test['validation_criteria']
		
		res = VerifySimResults( results_Neuron, results_Eden, criteria,
			skip_missing_criteria=skip_missing, explain_verification_on_fail=explain_verification_on_fail)
		if res is True:
			ok = True
			print('PASS %s %s' % ( test_filename, testcase_type ) )
	else:
		raise ValueError('Unknown test type "'+ typ + '"')
	
	return ok

def RunTests(
	tests,
	skip_missing_criteria=False, explain_verification_on_fail=False,
	extra_cmdline_args=None,
	verbose = False,
):
	failed_tests = []
	for test in tests:
		ok = RunTest( test, skip_missing_criteria, explain_verification_on_fail, extra_cmdline_args, verbose )
		
		if ok is not True:
			failed_tests.append(test)

	if not failed_tests:
		return True
	else:
		return failed_tests
