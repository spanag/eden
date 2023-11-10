import sys
from eden_tools import *
from eden_simulator import runEden
from eden_tools.run_sim import runNeuron

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

test_nml_dir = 'neuroml/'
tests = [
{
	'type': 'smoke_test',
	'sim_file': test_nml_dir + 'LEMS_EdenTest_PassiveCompartment.xml',
},
{
	'type': 'half_vs_half',
	'sim_file': 'neuroml/' + 'EdenTest_Extension_CustomSetup.nml',
	'validation_criteria':[
		["popA/0/izTonicSpiking/v"  , "pop0[0]/v", "exact"],
		["popB/0/izPhasicSpiking/v" , "pop0[1]/v", "exact"],
		["popC/0/izTonicBursting/v" , "pop0[2]/v", "exact"],
		["popD/0/izPhasicBursting/v", "pop0[3]/v", "exact"],
		["popE/0/izMixedMode/v"     , "pop0[4]/v", "exact"],
		["popp[0]/0/v"              , "popp[2]/0/v", "exact"],
		["popp[1]/0/v"              , "popp[3]/0/v", "exact"],
	]
},
{
	'type': 'half_vs_half',
	'sim_file': 'neuroml/' + 'EdenTest_Extension_VariableRequirement.nml',
	'validation_criteria':[
		["popA[0]/v", "pRp[0]/v", { 'type': 'box', 'dt': 0, 'dv': 0.0000005 }], # due to numerical roundoff really
		["popA[0]/u", "pSp[0]/u", { 'type': 'box', 'dt': 0, 'dv': 1.5e-18 }], # due to numerical roundoff really
		["popB[0]/v", "pRp[1]/v", "exact"],
		["popB[0]/u", "pSp[1]/u", { 'type': 'box', 'dt': 0, 'dv': 1.5e-18 }], # due to numerical roundoff really
	]
},
{
	'type': 'half_vs_half',
	'sim_file': 'neuroml/' + 'EdenTest_Extension_CustomIO.nml',
	'validation_criteria':[
		["pop0[0]/v", "popX[0]/v", "exact"],
		["pop0[1]/v", "popX[1]/v", "exact"],
		["pop0[2]/v", "popX[2]/v", { 'type': 'box', 'dt': 0, 'dv': 0.000006 }], # due to numerical roundoff really
        ["pop0[3]/v", "popX[3]/v", "exact"],
	]
},
{
	'sim_file': 'neuroml/' + 'EdenTest_Extension_WritableRequirement.nml',
	'type': 'half_vs_half',
	'validation_criteria': [
		["proj0[0]/post/flag", "pop0[0]/flag", "exact"],
	],
},
{
	'type': 'validate_vs_neuron',
	'sim_file': test_nml_dir + 'LEMS_EdenTest1.xml',
	'validation_criteria': {
		'pop0/0/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0005, 'dv': 0.002 },
		'pop0/1/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0005, 'dv': 0.0001 },
		'pop0/2/MultiCompCell/1/v': { 'type': 'box', 'dt': 0.0005, 'dv': 0.002 },
		'pop0/2/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0005, 'dv': 0.002 },
		'pop0/2/MultiCompCell/1/bioPhys1/membraneProperties/naChans/naChan/m/q': { 'type': 'box', 'dt': 0.0005, 'dv': 0.005 },
		'pop0/2/MultiCompCell/0': {'type': 'samecount_box', 'dt': 0.0002 },
	},
},
{
	'type': 'smoke_test',
	'sim_file': test_nml_dir + 'LEMS_EdenTest2.xml',
},
{
	'type': 'smoke_test',
	'sim_file': test_nml_dir + 'LEMS_EdenTest3.xml',
},
{
	'type': 'smoke_test',
	'sim_file': test_nml_dir + 'LEMS_EdenTest4.xml',
},
{
	'type': 'validate_vs_neuron',
	'sim_file': test_nml_dir + 'LEMS_EdenTest_CoreInputs.xml',
	'validation_criteria': {
		'InputDemos/01/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000010, 'dv': 0.000100 },
		'InputDemos/02/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000010, 'dv': 0.001000 },
		# 'InputDemos/03/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000001, 'dv': 0.000100 },
		# 'InputDemos/04/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000001, 'dv': 0.000100 },
		# 'InputDemos/05/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000001, 'dv': 0.000100 },
		'InputDemos/06/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000020, 'dv': 0.000500 },
		# 'InputDemos/07/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000001, 'dv': 0.000100 },
		# 'InputDemos/08/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000001, 'dv': 0.000100 },
		'InputDemos/09/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000010, 'dv': 0.000100 },		
		'InputDemos/11/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000010, 'dv': 0.000100 },
		'InputDemos/12/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000500, 'dv': 0.000500 },
		# 'InputDemos/13/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000001, 'dv': 0.000100 },
		# 'InputDemos/14/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000001, 'dv': 0.000100 },
		# 'InputDemos/15/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000001, 'dv': 0.000100 },
		# 'InputDemos/16/PassiveCell/0/v': { 'type': 'box', 'dt': 0.000001, 'dv': 0.000100 },
		# 'InputDemosDL[01]/v': { 'type': 'box', 'dt': 0.000001, 'dv': 0.000100 },
		# 'InputDemosDL[02]/v': { 'type': 'box', 'dt': 0.000001, 'dv': 0.000100 },
		# 'InputDemosDL[03]/v': { 'type': 'box', 'dt': 0.000001, 'dv': 0.000100 },
		
		'Input_spikeArray[0]'    : {'type': 'samecount_box', 'dt': 0.00001001 },
		'Input_spikeGenerator[0]': {'type': 'samecount_box', 'dt': 0.00001001 },
		
		# TODO automate running the ones commented out with jLEMS
	},
},
{
	'type': 'validate_vs_neuron',
	'sim_file': test_nml_dir + 'LEMS_EdenTest_CoreSynapses.xml',
	'validation_criteria': {
		'ChemPre/0/HHCell/0/v'              : { 'type': 'box', 'dt': 0.000050, 'dv': 0.000200 },

		'ChemPost/01/PassiveCell/0/v'       : { 'type': 'box', 'dt': 0.000010, 'dv': 0.000100 },
		'ChemPost/02/PassiveCell/0/v'       : { 'type': 'box', 'dt': 0.000020, 'dv': 0.000100 },
		'ChemPost/03/PassiveCell/0/v'       : { 'type': 'box', 'dt': 0.000020, 'dv': 0.000200 },
		'ChemPost/04/PassiveCell/0/v'       : { 'type': 'box', 'dt': 0.000020, 'dv': 0.000200 },
		'ChemPost/05/PassiveCell/0/v'       : { 'type': 'box', 'dt': 0.000020, 'dv': 0.000100 },
		'ChemPost/06/PassiveCell/0/v'       : { 'type': 'box', 'dt': 0.000010, 'dv': 0.000100 },
		'ChemPost/07/PassiveCell/0/v'       : { 'type': 'box', 'dt': 0.000030, 'dv': 0.000100 },
		'ChemPost/08/PassiveCell/0/v'       : { 'type': 'box', 'dt': 0.000030, 'dv': 0.000100 },
		'ChemPost/09/PassiveCell/0/v'       : { 'type': 'box', 'dt': 0.000030, 'dv': 0.000100 },
		'ChemPost/10/PassiveCell/0/v'       : { 'type': 'box', 'dt': 0.000020, 'dv': 0.000500 },
		
		'ChemPost/21/PassiveCell/0/v'       : { 'type': 'box', 'dt': 0.000010, 'dv': 0.000100 },
		
		
		'ContPre/01/PassiveCell/0/v'        : { 'type': 'box', 'dt': 0.000010, 'dv': 0.000100 },
		'ContPre/02/PassiveCell/0/v'        : { 'type': 'box', 'dt': 0.000010, 'dv': 0.000100 },
		'ContPre/03/PassiveCell/0/v'        : { 'type': 'box', 'dt': 0.000010, 'dv': 0.000100 },
		
		'ContPost/02/PassiveCell/0/v'       : { 'type': 'box', 'dt': 0.000010, 'dv': 0.000100 },
		'ContPost/03/PassiveCell/0/v'       : { 'type': 'box', 'dt': 0.000010, 'dv': 0.000100 },
		'ContPost/04/PassiveCell/0/v'       : { 'type': 'box', 'dt': 0.000010, 'dv': 0.000100 },
	},
},
{
	'type': 'validate_vs_neuron',
	'sim_file': test_nml_dir + 'LEMS_EdenTest_ArtificialCells.xml',
	'validation_criteria': {
		'Pop_iafCell[0]/v'       : { 'type': 'box', 'dt': 0.000100, 'dv': 0.000100 },
		'Pop_iafTauCell[0]/v'    : { 'type': 'box', 'dt': 0.000100, 'dv': 0.000100 },
		'Pop_iafRefCell[0]/v'    : { 'type': 'box', 'dt': 0.000200, 'dv': 0.000100 },
		'Pop_iafTauRefCell[0]/v' : { 'type': 'box', 'dt': 0.000100, 'dv': 0.000100 },
		# 'Pop_izhikevichCell[0]/v': { 'type': 'box', 'dt': 0.000500, 'dv': 0.005000 },
		'Pop_izhikevich2007Cell[0]/v': { 'type': 'box', 'dt': 0.000100, 'dv': 0.000500 },
		# 'Pop_adExIaFCell[0]/v' : { 'type': 'box', 'dt': 0.000100, 'dv': 0.000100 },
		'Pop_fitzHughNagumoCell[0]/V':     { 'type': 'box', 'dt': 0.000001, 'dv': 0.000100 },
		'Pop_fitzHughNagumo1969Cell[0]/V': { 'type': 'box', 'dt': 0.000001, 'dv': 0.000100 },
		# 'Pop_pinskyRinzelCA3Cell[0]/Vs': { 'type': 'box', 'dt': 0.000020, 'dv': 0.000100 },
		'Pop_IF_curr_alpha[0]/v' : { 'type': 'box', 'dt': 0.000100, 'dv': 0.000100 },
		'Pop_IF_curr_exp[0]/v'   : { 'type': 'box', 'dt': 0.000100, 'dv': 0.000100 },
		'Pop_IF_cond_alpha[0]/v' : { 'type': 'box', 'dt': 0.000100, 'dv': 0.000100 },
		'Pop_IF_cond_exp[0]/v'   : { 'type': 'box', 'dt': 0.000100, 'dv': 0.000100 },
		# 'Pop_EIF_cond_exp_isfa_ista[0]/v'   : { 'type': 'box', 'dt': 0.000100, 'dv': 0.000800 },
		# 'Pop_EIF_cond_alpha_isfa_ista[0]/v' : { 'type': 'box', 'dt': 0.000100, 'dv': 0.000100 },
		'Pop_HH_cond_exp[0]/v'   : { 'type': 'box', 'dt': 0.000020, 'dv': 0.000100 },
		# 'Pop_iafTauCell_Quadratic[0]/v'     : { 'type': 'box', 'dt': 0.000020, 'dv': 0.000100 },
		'Pop_iafTauCell_Quadratic[0]/v'     : { 'type': 'box', 'dt': 0.000020, 'dv': 0.000100 },
		
		"Pop_iafTauCell[0]"        :  {'type': 'samecount_box', 'dt': 0.00001001 },
		"Pop_iafTauRefCell[0]"     :  {'type': 'samecount_box', 'dt': 0.00002001 },
		"Pop_iafCell[0]"           :  {'type': 'samecount_box', 'dt': 0.00001001 },
		"Pop_iafRefCell[0]"        :  {'type': 'samecount_box', 'dt': 0.00010001 },
		# "Pop_izhikevichCell[0]"    :  {'type': 'samecount_box', 'dt': 0.00001001 },
		"Pop_izhikevich2007Cell[0]":  {'type': 'samecount_box', 'dt': 0.0015001 },
		# "Pop_adExIaFCell[0]"       :  {'type': 'samecount_box', 'dt': 0.00001001 },
		"Pop_IF_curr_alpha[0]"     :  {'type': 'samecount_box', 'dt': 0.00001001 },
		"Pop_IF_curr_exp[0]"       :  {'type': 'samecount_box', 'dt': 0.00001001 },
		"Pop_IF_cond_alpha[0]"     :  {'type': 'samecount_box', 'dt': 0.00001001 },
		"Pop_IF_cond_exp[0]"       :  {'type': 'samecount_box', 'dt': 0.00001001 },
	},
},
{
	'type': 'validate_vs_neuron',
	'sim_file': test_nml_dir + 'LEMS_EdenTest_CableDiscretisation.xml',
	'validation_criteria': {
		'Pop_Unified/0/Cable_Unified/0/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.003000 },
		'Pop_Unified/0/Cable_Unified/1/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.003000 },
		'Pop_Unified/0/Cable_Unified/2/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.003000 },
		'Pop_Subdivided/0/Cable_Subdivided/0/v'  : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
		'Pop_Subdivided_Explicitly/0/Cable_Subdivided_Explicitly/0/v'      : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
		'Pop_Subdivided_Explicitly/0/Cable_Subdivided_Explicitly/1/v'      : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
		'Pop_Subdivided_Explicitly/0/Cable_Subdivided_Explicitly/2/v'      : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
		'Pop_Subdivided_Explicitly/0/Cable_Subdivided_Explicitly/3/v'      : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
		'Pop_Subdivided_Explicitly/0/Cable_Subdivided_Explicitly/4/v'      : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
		'Pop_Subdivided_Effectively/0/Cable_Subdivided_Effectively/0/v'    : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
		'Pop_Subdivided_Effectively/0/Cable_Subdivided_Effectively/1/v'    : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
		'Pop_Subdivided_Effectively/0/Cable_Subdivided_Effectively/2/v'    : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
		'Pop_Subdivided_Effectively/0/Cable_Subdivided_Effectively/3/v'    : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
		'Pop_Subdivided_Effectively/0/Cable_Subdivided_Effectively/4/v'    : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
	},
},
{
	'type': 'validate_vs_neuron',
	'sim_file': test_nml_dir + 'LEMS_EdenTest_InhomogeneousParameters.xml',
	'validation_criteria': {
		'Pop_Flat/0/SimpleCell_Flat/0/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
		'Pop_Flat/0/SimpleCell_Flat/1/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.002000 },
		'Pop_Flat/0/SimpleCell_Flat/2/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.002000 },
		'Pop_InhF/0/SimpleCell_InhF/0/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
		'Pop_InhF/0/SimpleCell_InhF/1/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.002000 },
		'Pop_InhF/0/SimpleCell_InhF/2/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.002000 },
		'Pop_Inho/0/SimpleCell_Inho/0/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
		'Pop_Inho/0/SimpleCell_Inho/1/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.002000 },
		'Pop_Inho/0/SimpleCell_Inho/2/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.002000 },
		'Pop_InhN/0/SimpleCell_InhN/0/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.001000 },
		'Pop_InhN/0/SimpleCell_InhN/1/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.002000 },
		'Pop_InhN/0/SimpleCell_InhN/2/v'        : { 'type': 'box', 'dt': 0.001000, 'dv': 0.002000 },
	},
},
{
	'type': 'validate_vs_neuron',
	'sim_file': test_nml_dir + 'LEMS_EdenTest_DomainDecomposition.xml',
	'validation_criteria': {
		'pop0/00/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0015, 'dv': 0.004 },
		'pop0/01/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0015, 'dv': 0.004 },
		'pop0/02/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0015, 'dv': 0.004 },
		'pop0/03/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0015, 'dv': 0.004 },
		'pop0/04/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0015, 'dv': 0.004 },
		'pop0/05/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0015, 'dv': 0.004 },
		'pop0/06/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0015, 'dv': 0.004 },
		'pop0/07/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0015, 'dv': 0.004 },
		'pop0/08/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0015, 'dv': 0.004 },
		'pop0/09/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0015, 'dv': 0.004 },
		'pop0/10/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0015, 'dv': 0.004 },
		'pop0/11/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0015, 'dv': 0.004 },
		'pop0/12/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0015, 'dv': 0.004 },
		'pop0/13/MultiCompCell/0/v': { 'type': 'box', 'dt': 0.0015, 'dv': 0.004 },
	},
},
{
	'type': 'eden_vs_eden',
	'sim_file': test_nml_dir + 'LEMS_EdenTest_DomainDecomposition.xml',
	# 'truth_kwargs': {},
	'test_kwargs': { 'full_cmdline': ['mpirun','-n','4','eden-mpi', 'nml', test_nml_dir + 'LEMS_EdenTest_DomainDecomposition.xml' ], 'threads':2, 'verbose': True },
	'validation_criteria': 'exact'
},

]

# add some auxiliary processes, just to test pipes for now
import subprocess
concurrent_processes = []
fifo_pair_filenames = [('event', 'spikeDelayed.from_sim.pipe','spikeDelayed.to_sim.pipe'), ('trajectory', 'timeseriesDelayed.from_sim.pipe','timeseriesDelayed.to_sim.pipe'),]
runtime_dir = '.'
try:
	for type, s_a, a_s in fifo_pair_filenames:
		s_a, a_s = ( f'{runtime_dir}/'+x for x in (s_a, a_s))
		subprocess.call(['rm', s_a, a_s]);
		subprocess.call(['mkfifo', s_a, a_s]);
		p = subprocess.Popen(['python', 'delayeur.py', '0.010', type, f'{s_a}', f'{a_s}'])
		concurrent_processes.append(p)
	
	# now run the tests
	res = RunTests(tests, verbose = True)
	if res is True:
		print('All tests passed')
		sys.exit(0)
	else:
		print('%d tests failed !' % len(res))
		sys.exit(1)
	
finally:
	# clean up after the processes
	for p in concurrent_processes:
		try:
			p.wait(timeout=0.1) # give some time for graceful shutdown of aux processes
		except subprocess.TimeoutExpired: pass
	if any(p for p in concurrent_processes if p.returncode is None): print('Terminating remaining processes...')
	for p in (p for p in concurrent_processes if p.returncode is None):
		print(p)
		p.kill()
