import sys
from eden_tools import *

test_nml_dir = 'neuroml/'
tests = [
{
	'type': 'smoke_test',
	'sim_file': test_nml_dir + 'LEMS_EdenTest_PassiveCompartment.xml',
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
		
		# TODO automate running the ones commented out with jLEMSa
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
res = RunTests(tests, verbose = True)
if res is True:
	print('All tests passed')
	sys.exit(0)
else:
	print('%d tests failed !' % len(res))
	sys.exit(1)
