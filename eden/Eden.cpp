/*
       ###########    ###############              ###########    ######     #########      
     ###############   ##################        ###############    #####  ##############   
   #####          ####      #####    ######    #####          ####     ######       ######  
  ####                    ####          ####  ####                     ####           ##### 
  ###############         ###            ###  ###############          ####            #### 
  ############           ####            ###  ############             ###             #### 
  ####                   ####            ###  ####                     ###             #### 
  ####              #    ####           ####  ####              #     ####             #### 
   #####          ###     ####         ####    #####          ###     ####            ####  
    ################        ##############      ################     ####            ####   
      ##########              ##########          ##########        ####            #####   
*/

/*
Extensible Dynamics Engine for Networks
Parallel simulation engine for ODE-based models
*/

#include "Common.h"
#include "NeuroML.h"

#include <math.h>
#include <limits.h>
#include <errno.h>


#include <map>
#include <set>
#include <chrono>

#include "MMMallocator.h"

#include <omp.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#if defined (__linux__) || defined(__APPLE__)
#include <dlfcn.h> // for dynamic loading
#endif

#ifdef _WIN32
// for dynamic loading and other OS specific stuff
// #include <windows.h> // loaded through Common.h at the moment, TODO break out in Windows specific header 
#endif

extern "C" {	

// do not specify alignment for the pointers, in the generic interface
// cannot specify __restrict__ because it is quietly dropped by compilers ( ! ) when the type is allocated with new
// causing a type mismatch when operator delete(T * __restrict__) is called (then why didn't they drop __restrict__ from there too ??)
// just hope the type "mismatch" won't cause a crash in practice
typedef float * Table_F32;
typedef long long * Table_I64;

// assume standard C calling convention, which is probably the only cross-module one in most architectures
// bother if problems arise LATER
typedef void ( *IterationCallback)(
	double time,
	float dt,
	const float *__restrict__ constants, long long const_local_index,
	const long long *__restrict__ const_table_f32_sizes, const Table_F32 *__restrict__ const_table_f32_arrays, long long table_cf32_local_index,
	const long long *__restrict__ const_table_i64_sizes, const Table_I64 *__restrict__ const_table_i64_arrays, long long table_ci64_local_index,
	const long long *__restrict__ state_table_f32_sizes, const Table_F32 *__restrict__ state_table_f32_arrays, Table_F32 *__restrict__ stateNext_table_f32_arrays, long long table_sf32_local_index,
	const long long *__restrict__ state_table_i64_sizes,       Table_I64 *__restrict__ state_table_i64_arrays, Table_I64 *__restrict__ stateNext_table_i64_arrays, long long table_si64_local_index,
	const float *__restrict__ state, float *__restrict__ stateNext, long long state_local_index,
	long long step
	
);
void UndefinedCallback(
	double time,
	float dt,
	const float *__restrict__ constants, long long const_local_index,
	const long long *__restrict__ const_table_f32_sizes, const Table_F32 *__restrict__ const_table_f32_arrays, long long table_cf32_local_index,
	const long long *__restrict__ const_table_i64_sizes, const Table_I64 *__restrict__ const_table_i64_arrays, long long table_ci64_local_index,
	const long long *__restrict__ state_table_f32_sizes, const Table_F32 *__restrict__ state_table_f32_arrays, Table_F32 *__restrict__ stateNext_table_f32_arrays, long long table_sf32_local_index,
	const long long *__restrict__ state_table_i64_sizes,       Table_I64 *__restrict__ state_table_i64_arrays, Table_I64 *__restrict__ stateNext_table_i64_arrays, long long table_si64_local_index,
	const float *__restrict__ state, float *__restrict__ stateNext, long long state_local_index,
	long long step
){
	// what else? LATER
	abort();
}
}


// helpers

template< typename Container >
static void AppendToVector(Container &append_to, const Container &append_this){
	append_to.insert(append_to.end(), append_this.begin(), append_this.end());
};
template<
	typename CAppendTo, typename CAppendThis,
	typename std::enable_if< !std::is_same< CAppendTo, CAppendThis >::value, int >::type = 0
>
static void AppendToVector(CAppendTo &append_to, const CAppendThis &append_this){
	auto new_size = append_to.size() + append_this.size() ;
	append_to.reserve( new_size );
	for( size_t i = 0; i < append_this.size(); i++ ){
		append_to.push_back( append_this[i] );
	}
};

// LATER add a helper for single variables

// do not use default accuracy when converting numerics to alpha !
// expicitly specify what the alpha is used for
template< typename T, typename = typename std::enable_if< std::is_integral<T>::value >::type >
std::string accurate_string( T val ){
	return std::to_string( val );
}
std::string accurate_string( float val ){
	char tmps[100];
	sprintf(tmps, "%.9g", val);
	return tmps;
}
std::string accurate_string( double val ){
	char tmps[100];
	sprintf(tmps, "%.17g", val);
	return tmps;
}

std::string presentable_string( double val ){
	char tmps[100];
	sprintf(tmps, "%g", val);
	return tmps;
}
template< typename T, typename = typename std::enable_if< std::is_integral<T>::value >::type >
std::string presentable_string( T val ){
	char tmps[100];
	if( val < -1000000 || 1000000 < val ){
		// extreme values are probably packed refs
		sprintf(tmps, "0x%llx", (long long) val);
	}
	else{
		sprintf(tmps, "%lld", (long long) val);
	}
	
	return tmps;
}

template< typename T, typename std::enable_if< std::is_integral<T>::value, int >::type = 0 >
std::string itos(T val){
	return accurate_string(val);
}
template< typename T, typename std::enable_if< std::is_enum<T>::value, int >::type = 0 >
std::string itos(T val){
	return std::to_string(val);
}

// for strictly fixed-width data output (also useful for parallel writing to files)
struct FixedWidthNumberPrinter{
	int column_size;
	int delimiter_size;
	char delimiter_char;
	char format[50];
	
	int getNumberSize() const {
		return column_size - delimiter_size;// for separator
	}
	
	FixedWidthNumberPrinter( int _csize, char _delcha = ' ', int _dellen = 1 ){
		column_size = _csize;
		delimiter_size = _dellen;
		delimiter_char = _delcha;
		assert( column_size > delimiter_size );
		
		const int number_size = getNumberSize();
		const int digits = column_size - 3 - 5; // "+1.", "e+308"
		sprintf(format, "%%+%d.%dg", number_size, digits );
	}
	
	// should have length of column_size + terminator
	void write( float val, char *buf) const {
		const int number_size = getNumberSize();
		snprintf( buf, number_size + 1, format, val );
		// Also add some spaces, to delimit columns
		for( int i = 0; i < delimiter_size; i++ ){
			buf[number_size + i] = delimiter_char;
		}
		buf[number_size+1] = '\0';
	}
};

//Linear interpolation from a -> b with respective parameter from 0 -> 1
static Real Lerp(Real a, Real b, Real ratio){
	return a + (b-a)*ratio;
}

template <typename Real>
struct GeomHelp_Base{
	inline static Real Length(Real dx, Real dy, Real dz){
		return std::sqrt(dx*dx + dy*dy + dz*dz);
	}
	inline static Real Area(Real length, Real diam_proximal, Real diam_distal){
		// actually just the external surface area of the frustum
		// what about the end butt of dendrites (and soma which is actually large!), though?
		// do what NEURON does, ignore them and let higher-level modelling software apply corrections
		if(length == 0){
			// spherical soma or something, TODO more robust detection at parse time too!
			// XXX this is ambiguous, since a section with zero length can be thought of as an annular disc!!
			// A workaround for the current use is to add a minuscule non-zero length to segments that represent annular discs, to follow the appropriate code path  for that case
			return M_PI * diam_distal * diam_distal;
		}
		else return (M_PI / 2.0) * (diam_proximal + diam_distal) * std::sqrt( ((diam_proximal - diam_distal) * (diam_proximal - diam_distal) / 4.0) + length * length);
	}
	inline static Real Volume(Real length, Real diam_proximal, Real diam_distal){
		if(length == 0){
			// spherical soma or something, TODO more robust detection at parse time too!
			return (M_PI / 6.0) * diam_distal * diam_distal * diam_distal;
		}
		else return (M_PI / 3.0) * length * ( diam_proximal*diam_proximal + diam_distal*diam_distal + diam_proximal*diam_distal ) / 4.0;
	}
};
typedef GeomHelp_Base<float> GeomHelp;

// some math on (x,y,z,d) tuples of segments, to avoid copy pasting 
static Morphology::Segment::Point3DWithDiam LerpPt3d(
	Morphology::Segment::Point3DWithDiam from,
	Morphology::Segment::Point3DWithDiam to, Real ratio
){
	Morphology::Segment::Point3DWithDiam ret = {
		Lerp(from.x, to.x, ratio),
		Lerp(from.y, to.y, ratio),
		Lerp(from.z, to.z, ratio),
		Lerp(from.d, to.d, ratio)
	};
	return ret;
}
static Real DistancePt3d(
	Morphology::Segment::Point3DWithDiam from,
	Morphology::Segment::Point3DWithDiam to
){
	return GeomHelp::Length( to.x - from.x, to.y - from.y, to.z - from.z );
}

// MPI context, just a few globals like world size, rank etc.
#ifdef USE_MPI
struct MpiContext{
	int world_size;
	int rank;
};
const int MYMPI_TAG_BUF_SEND = 99;
MpiContext my_mpi;


// TODO use logging machinery in whole codebase
FILE *fLog = stdout;
void Say( const char *format, ... ){
	va_list args;
	va_start(args, format);
	
	std::string new_format = "rank "+std::to_string(my_mpi.rank)+" : " + format + "\n";
	vfprintf(fLog, new_format.c_str(), args);
	fflush(stdout);
	
	va_end (args);
}

#endif


// initial states, internal constants, connectivity matrices, iteration function pointers and everything
// so crunching can commence
struct RawTables{
	
	const static size_t ALIGNMENT = 32;
	
	typedef std::vector< float, _mm_Mallocator<float, ALIGNMENT> > Table_F32;
	typedef std::vector< long long, _mm_Mallocator<long long, ALIGNMENT> > Table_I64;
	
	// TODO aligned vectors, e.g. std::vector<T, boost::alignment::aligned_allocator<T, 16>>
	Table_F32 global_initial_state; // sum of states of work items
	Table_F32 global_constants; // sum of states of work items
	Table_I64 index_constants; // TODO if ever
	
	std::vector<long long> global_state_f32_index; // for each work unit TODO explain they are not completely "global"
	std::vector<long long> global_const_f32_index; // for each work unit
	
	//the tables TODO aligned
	std::vector<long long> global_table_const_f32_index; // for each work unit
	std::vector<long long> global_table_const_i64_index; // for each work unit
	std::vector<long long> global_table_state_f32_index; // for each work unit
	std::vector<long long> global_table_state_i64_index; // for each work unit
	//using std::vector due to a-priori-unknown size of tables; will use a pointer-based allocator to pass to callbacks as vectors, and compaction LATER
	std::vector<Table_F32> global_tables_const_f32_arrays; //the backing store for each table
	std::vector<Table_I64> global_tables_const_i64_arrays; //the backing store for each table
	std::vector<Table_F32> global_tables_state_f32_arrays; //the backing store for each table
	std::vector<Table_I64> global_tables_state_i64_arrays; //the backing store for each table
	
	std::vector<IterationCallback> callbacks; // for each work unit
	
	// some special-purpose tables
	
	// These are to access the singular, flat state & const vectors. They are not filled in otherwise.
	long long global_const_tabref;
	long long global_state_tabref;
	
	// the data loggers use this to access the variables, for remapping to be possible in MPI case - why not in eng.config then?
	// long long logger_table_const_i64_to_state_f32_index;
	
	// TODO what about the MPI mirrors ?
	
	RawTables(){
		global_const_tabref = -1;
		global_state_tabref = -1;
	}
};
// references to the raw tables
struct TabEntryRef{
	long long table;
	int entry;
	// could also use the compressed, encoded combination
};
typedef long long TabEntryRef_Packed;
auto GetEncodedTableEntryId = []( long long global_idx_T_dest_table, long long entry_idx_T_dest ){
	// pack 1 trillion tables -> 16 million entries into 64bit indexes, upgrade if needed LATER
	
	const unsigned long long table_id = global_idx_T_dest_table * (1 << 24);
	const unsigned long long entry_id = entry_idx_T_dest % (1 << 24);
	const unsigned long long packed_id = table_id | entry_id ;
	
	return packed_id;
};
auto GetDecodedTableEntryId = [  ]( unsigned long long packed_id ){
	long long global_idx_T_dest_table = packed_id >> 24;
	int entry_idx_T_dest = packed_id % (1 << 24);
	TabEntryRef ret = { global_idx_T_dest_table, entry_idx_T_dest };
	return ret;
};

// smuggle int32 into f32 position
// dodgy, but that's life
union TypePun_I32F32{
	int32_t i32; float f32;
	static_assert( sizeof(i32) == sizeof(f32), "Single-precision float must have same same size as int32_t for type punning to work" );
};
auto EncodeI32ToF32( int32_t i ){
	TypePun_I32F32 cast;
	cast.i32 = i;
	return cast.f32;
}
auto EncodeF32ToI32( float f ){
	TypePun_I32F32 cast;
	cast.f32 = f;
	return cast.i32;
}


// and more information that is needed for the engine
struct EngineConfig{
	
	struct TrajectoryLogger {
		
		struct LogColumn{
			enum Type{
				NONE,
				TOPLEVEL_STATE,
				TABLE_STATE
			};
			enum ValueType{
				F32,
				I64
			};
			
			Type type;
			ValueType value_type;
			
			size_t entry;
			// size_t table;
			double scaleFactor; // in case of float values
			
			#ifdef USE_MPI
			int on_node;
			LogColumn(){
				on_node = -1;
			}
			#endif
		};
		
		std::string logfile_path;
		std::vector<LogColumn> columns;
		
	};
	
	
	long long work_items;
	double t_initial; // in engine time units
	double t_final;
	float dt; // in engine time units
	
	std::vector<TrajectoryLogger> trajectory_loggers;
	
	// for inter-node communication
	struct SendList_Impl{
		
		// std::vector<int> vpeer_positions_in_buffer;
		std::vector<size_t> vpeer_positions_in_globstate; // TODO packed table entries
		
		std::vector< TrajectoryLogger::LogColumn > daw_columns;
		// std::vector<int> daw_positions_in_buffer;
		// std::vector<size_t> daw_positions_in_globstate; // TODO packed table entries
		
		size_t spike_mirror_buffer;
	};
	struct RecvList_Impl{
		size_t value_mirror_buffer;
		ptrdiff_t value_mirror_size;
		// refs to mirror_buffer in table accesses, as well as off-table trajectory loggers, are resolved
		
		std::vector< std::vector<TabEntryRef_Packed> > spike_destinations;
	};
	std::map< int, SendList_Impl > sendlist_impls;
	std::map< int, RecvList_Impl > recvlist_impls;
	
	//spikes are triggered in buffer, they're automatically gathered
	
};


// general options for the simulator
struct SimulatorConfig{
	
	bool override_random_seed;
	long long override_random_seed_value;
	
	bool verbose;// TODO implement in more places	
	
	bool debug;
	bool debug_netcode;
	
	bool dump_raw_state_scalar;
	bool dump_raw_state_table;
	bool dump_raw_layout;
	
	bool use_icc;
	bool tweak_lmvec;
	bool output_assembly;
	
	// TODO knobs:
	// vector vs.hardcoded sequence for bwd euler, also heuristic
	enum CableEquationSolver{
		CABLE_SOLVER_AUTO,
		CABLE_FWD_EULER,
		CABLE_BWD_EULER,
	};
	CableEquationSolver cable_solver;
	
	SimulatorConfig(){
		verbose = false;
		debug = false;
		debug_netcode = false;
		
		dump_raw_state_scalar = false;
		dump_raw_state_table = false;
		dump_raw_layout = false;
		
		use_icc = false;
		tweak_lmvec = false;
		output_assembly = false;
		
		cable_solver = CABLE_SOLVER_AUTO;
		
		override_random_seed = false;
	}
};

bool GenerateModel(const Model &model, const SimulatorConfig &config, EngineConfig &engine_config, RawTables &tabs){
	
	/*
	TODO:
		decouple derivative from integration rule
		split large cells into compartment-sized work units
	*/
	
	//------------------>  GENERAL INFORMATION
	/*
	Simulation parallelism assumes calculations are split in essential units, and calculations for each unit are themselves computed sequentially.
	These units may be, for example, discretized elements in a PDE model, or functional units that interact with each other.
	Each unit should include all calculations needed for a state variable, to avoid splitting the calculations in 
	multiple interlocked stages, which would add synchronization overhead.
	In some architectures, introducing a flux-gathering stage before the state-progression stage might perhaps be called for, 
	due to the irregular flux connectivity between units (each unit may have a different degree of connectivity, which could complicate things otherwise)
	
	Sometimes, a work unit may involve a set of constants and/or state variables, whose size varies across instances of the same work unit type.
	In particular, the size may be different for each work unit, so the work unit cannot be specialized into a fixed number of working set sizes.
	This is the case for example, for irregular graph models where the work unit is a single vertex.
	In this case, a distinction is made between the working set that is common for all work units, and the variable-sized part each instance has.
	The common working set is named the internal working set, and the variable-sized part is split into tables:
		Each table is defined by a data type, length and starting memory position.
		A variable-sized part of a work unit may work with one or more tables. A work unit may have multiple variable-sized parts.
		The set of tables to be used is the same for all work units of the same type.
			(If tables are empty or fixed-length or fixed-structure in most instances, the work unit could be specialized in sub-cases)
	
	
	The internal working set of all units is joined into single vectors, one for each value type (state, constant or index).
	The working set tables for each unit may be allocated anywhere, possibly joined along with the internal working sets, with each other, or anything really.
	
	// TODO explain discrete events.
	Each unit:
		writes its own set of state variables
		reads  its own set of state variables
		reads  its own set of parameters
		reads  global parameters
		reads  state variables of other units
		reads  global state variables
	using computation code unique to the type it is instantiated from.
	
	*/
	
	//------------------>  EDEN/NEUROML INFORMATION
	/*
	Internal states are located within compartments of neurons.
	State variables are:
		Compartment membrane voltages
		Intra/extracellular ion concentrations 
		Internal states of ion channels 
		Internal states of synapses
		Internal states of input components
	
	For signalling purposes, the firing state of each event port is also preserved as a state variable.
	
	Membrane potential changes at a rate of total inward current, divided by total compartment capacitance.
	Current enters(+) and leaves(-) a compartment through:
		channels, active or passive ( Ichan = gchan * V - Vextra)
		intra-cellular leaks between adjacent compartments (Iaxial)
		synapses, chemical (current moves to extracellular fluid) or electrical (current moves to adjacent cell)
		electrodes, artificial inputs? TODO specially handle voltage clamps!
	Ions enter and leave a compartment through:
		channels, active or passive
		artificial injection perhaps?
	A channel's gates or Markov states change at a rate influenced by:
		the channel's states
		compartment voltage
		certain ion or signaller concentrations
	A synapse's state may change at a rate influenced by:
		its state
		pre-synaptic membrane potential
		post-synaptic membrane potential
		spikes travelling through the synapses
	*/
	
	/*
	Assume the essential unit is a neuron, or a compartment.
	The factors differentiating calculations for each unit are:
		Internal mechanisms (ion channels, ion pools, selective longitudinal diffusion)
		Attached mechanisms (due to synapses and input sources)
	*/
	
	const auto &dimensions          = model.dimensions          ;
	const auto &component_types     = model.component_types     ;
	const auto &morphologies        = model.morphologies        ;
	const auto &biophysics          = model.biophysics          ;
	const auto &ion_species         = model.ion_species         ; (void)ion_species;
	const auto &conc_models         = model.conc_models         ;
	const auto &ion_channels        = model.ion_channels        ;
	const auto &cell_types          = model.cell_types          ;
	const auto &synaptic_components = model.synaptic_components ;
	const auto &input_sources       = model.input_sources       ;
	const auto &networks            = model.networks            ;
	const auto &simulations         = model.simulations         ;
	const auto &target_simulation   = model.target_simulation   ;
	
	const Simulation &sim = simulations.get(target_simulation);
	
	
	// the basic RNG seed.
	// Modify using mose sim properties, LATER
	long long simulation_random_seed;
	if( config.override_random_seed ){
		simulation_random_seed = config.override_random_seed_value;
	}
	else{
		if( sim.seed_defined ){
			simulation_random_seed = sim.seed;
		}
		else{
			// use seconds from 1970-01-01 00:00
			// TODO pick the random number somewhere more convenient
			auto time_now = std::chrono::system_clock::now();
			simulation_random_seed = std::chrono::duration_cast<std::chrono::seconds>(time_now.time_since_epoch()).count();
		}
	}
	
	
	const Network &net = networks.get(target_simulation);
	
	struct CellInternalSignature{
		// LATER perhaps optimize for tiny table sizes; for larger tables, the extra pointer, memory fragmentation, etc. are amortized costs
		// LATER find out if all tables can be consolidated into a single vector, to remove data structure complexity
		// see also:
		//  The design of a template structure for a generalized data structure definition facility, ICSE '76
		//	https://dl.acm.org/citation.cfm?id=807713
		
		struct TableInfo{
			std::string _description; // may be strignified from more important parameters, such as corresponding mechanisms LATER
			std::string Description() const {
				return _description;
			}
			TableInfo(const std::string _desc){
				_description = _desc;
			}
		};
		
		// ComponentLayout: A component type has a unique abstract layout (no. of constants and states).
		// ComponentValueInstance: A component instance may have a different set of constants and/or initial states.
		// ComponentSubSignature: An allocation of a single/multi component type consists of the mapping between abstract layout and subset of a work item signature values/tables respectively.
		
		// To realize a component instance, an allocation must have been made, and the value instance applied on the allocation.
		// To realize a single instance:
		//	the layout must be extracted
		// 	the value instance must be determined for the specific instance
		// 	the subsig must be allocated, its singular contents initialized with the value instance
		// To realize a single instance:
		//	the layout must be extracted
		// 	the subsig must be allocated once
		// 	the contents for each instance initialized with the value instance of its prototype
		// the above is yet TODO
		
		// a component may need to keep the values to clone them on the same allocation of subsignature
		struct ComponentValueInstance{
			std::vector<Real> properties;
			std::vector<Real> statevars;
		};
		// merge
		struct ComponentSubSignature{
			
			struct Entry{
				enum ValueType{
					UNSET,
					F32,
					I64
				};
				
				size_t index;
				ValueType type;
				
				Entry(){ type = UNSET; }
				Entry( size_t _i, int _t ) : index(_i), type((ValueType)_t) { } // ICC crashes if _t is a ValueType
			};
			
			//Int component;
			
			// Indexes that could either refer to the cell(or compartment)'s scalar signature , or new spawned tables
			// (depending on the nature of the component
			std::vector<Entry> properties_to_constants;
			std::vector<Entry> statevars_to_states;
			
		};
		
		
		// The mapping of properties to offsets, for common components
		struct SynapticComponentImplementation{
			// Spike-receiving synaptic components require:
			//	- a trigger state table for each comp.type, in post-synaptic cell
			//	- an entry referring to trigger state table, in pre-synaptic cell send spike table 
			// Vpeer-dependent synaptic components require:
			//	- a Vpeer table for each comp.type, with entries for voltage of each peer
			
			// constants
			
			size_t Table_Weight; // for all of them
			
			// optional, is -1 if not used
			size_t Table_Delay; // for spiking synapses & hybrids
			size_t Table_Vpeer; // for gap junctions    & hybrids
			
			// if a native type
				size_t Table_Erev;
				size_t Table_Gbase;
				size_t Table_Gbase2;
				
				size_t Table_Tau;
				size_t Table_Tau2;
				size_t Table_Tau3;
				
				size_t Table_Ibase;
			
			// and states
			
			size_t Table_Trig;  // for spiking synapses & hybrids
			size_t Table_NextSpike;  // pending spike queue for spiking synapses & hybrids. LATER will hold more than one spike.
			
			// if a native type
				size_t Table_Grel;
				size_t Table_Grel2;
			
			
			// constants and states, if a LEMS component
				ComponentSubSignature synapse_component;
				
			// more LEMS components if a blocking/plastic synapse
				ComponentSubSignature block_component;
				ComponentSubSignature plasticity_component;
			
			SynapticComponentImplementation(){
				Table_Weight = -1;
				Table_Delay = -1;
				
				Table_Gbase = -1;
				Table_Gbase2 = -1;
				Table_Tau = -1;
				Table_Tau2 = -1;
				Table_Tau3 = -1;
				Table_Ibase = -1;
				
				Table_Trig = -1;
				Table_NextSpike = -1;
				Table_Grel = -1;
				Table_Grel2 = -1;
			}
			
		};
		
		struct InputImplementation{
			
			// for all of them
			size_t Table_Weight;
			
			// if a native type, a few of these
			size_t Table_Imax;
			
			size_t Table_Duration;
			size_t Table_Delay; // could also be called 'start time'
			size_t Table_Period;
			
			size_t Table_Phase;
			size_t Table_Istart;
			size_t Table_Iend;
			
			// if a spike list
			size_t Table_SpikeListTimes; // merged between instances, split by NaN sentinels
			//size_t Table_SpikeListStarts; // constants, one for each instance
			size_t Table_SpikeListPos; // one for each instance, as they traverse their respective lists
			// could also add a "next spike time" cache to eliminate most random access TODO
			
			// if a firing synapse, just for its internals (not weight, trigger etc.)
			SynapticComponentImplementation synimpl;
			
			// if a LEMS component
			//ComponentValueSignature input_component_values;
			ComponentSubSignature component; 
			
			InputImplementation(){
				Table_Weight = -1;
			}
		};
		struct SpikeSendingImplementation{
			// if a native type
			ptrdiff_t Table_SpikeRecipients;
			// possibly delay in sender LATER, depends on formulation
			
			SpikeSendingImplementation(){
				Table_SpikeRecipients = -1;
			}
			// what TODO with multiple spike-sending things in a compartment ??
		};
		
		struct RngImplementation{
			ptrdiff_t Index_RngSeed; // to cf32 for now, move to integer LATER
			RngImplementation(){
				Index_RngSeed = -1;
			}
		};
		
		// The mapping of properties to offsets, for physical cells
		struct PhysicalCell{
			
			// The mapping of NeuroML segments to compartments.
			// When the the "cable" attribute is not used, each segment maps to one compartment.
			// But when cables are introduced, multiple segments could map to compartments and vice versa, and not exclusively:
			// 	two compartments could include (separate) parts of the same segment and two segments
			// This is because NeuroML2 <segment>s have the semantics of 3D points of a section, as well as sections themselves.
			struct CompartmentDiscretization{
				
				int32_t number_of_compartments;
				
				// for each segment(seg_seq), keep a list of the compartments it spans, and up to which fractionAlong it spans them.
				// NB: Each segment overlapping a cable is fully contained by it, thus it 
				std::vector< std::vector<int32_t> >segment_to_compartment_seq; 
				std::vector< std::vector<Real> >segment_to_compartment_fractionAlong;
				// for each compartment, keep a list of the segments it contains, and fractionAlong for beginning and end.
				// (Implicitly, each segment in between is included until its end and from its start.)
				std::vector< std::vector<int32_t> >compartment_to_segment_seq; 
				std::vector< Real >compartment_to_first_segment_start_fractionAlong;
				std::vector< Real >compartment_to_last_segment_end_fractionAlong;
				
				// Also keep the tree-parent relationships, for they are not immediately clear without re-analysing cables
				std::vector<int32_t>tree_parent_per_compartment; //-1 for root
				
				// Note that every cell could be unique, and perhaps 16 bits are enough to describe neurons...
				// LATER efficiency trick: keep the above lists empty, non existent or something, for cases where the mapping is 1 to 1
				
				
				
				// Get the comp_seq which contains the fractionAlong across seg_seq (for this abstract mapping; lauout in memory could be different)
				int32_t GetCompartmentForSegmentLocation(Int seg_seq, Real fractionAlong) const {
					// printf("Compartment for seg_seq %d fractionAlong %g \n", (int)seg_seq, fractionAlong );
					// There is a possibillity for fractionAlong to deviate due to roundoff, perhaps the value could be clipped instead
					assert( (0 <= fractionAlong) && (fractionAlong <= 1) );
					
					// the correct compartment is the first whose ending fractionAlong exceeds (thus it overlaps) the target fractionAlong
					// use std::lower_bound if you find it better, through measurement
					for( int i = 0; i < (int) segment_to_compartment_fractionAlong[seg_seq].size(); i++ ){
						// printf("\tcheck %d %g ...\n", (int)segment_to_compartment_seq[seg_seq][i] , segment_to_compartment_fractionAlong[seg_seq][i] );
						// <= used for fractionAlong = 1 case
						if( fractionAlong <= segment_to_compartment_fractionAlong[seg_seq][i] ){
							return segment_to_compartment_seq[seg_seq][i];
						}
					}
					// on second thought, the whole range of fractionAlong should be covered, since segment_to_compartment_fractionAlong[seg_seq] ends in 1...
					// Last chance fallback: produce the last segment (more likely to not include fractionAlong = 1, because of roundoff)
					assert(false);
					int32_t last_compartment = segment_to_compartment_fractionAlong[seg_seq].size() - 1;
					if(last_compartment >= 0){
						return last_compartment;
					}
					else return -1;
				}
				
			};
			CompartmentDiscretization compartment_discretization;
			
			// some indexing structure to map symbolic representation to realized memory addresses and vice versa
			// indices to retrieve which state variable is which
			// LATER replace ints by something more type-safe (such as ints for F32 tables)
			size_t Index_Voltages; // where are voltages for each segment (if compact cell) ? -1 for non-compact cell
			// TODO place everything in use, here (like areas and resistances) 
			
			enum CompartmentGrouping{
				AUTO,
				FLAT,
				GROUPED
			};
			CompartmentGrouping compartment_grouping;
			
			struct CompartmentGroupingImplementation{
				// definition
				std::vector< IdListRle > distinct_compartment_types; // comp_seq per compartment type 
				std::vector< std::string > preupdate_codes, postupdate_codes; // code blocks per compartment type
				
				
				// vectors per compartment, for now; could join per-list for less pointer chasing
				std::vector< int > r_off; // RNG offsets
				std::vector< int > c_off, s_off;
				std::vector< int > cf32_off, sf32_off, ci64_off, si64_off;
				
				// indices to grouped lists of compartments
				std::vector< size_t > Index_CompList; // per compartment type
				
				// indices to the allocated vectors
				size_t Index_Roff;
				size_t Index_Coff, Index_Soff;
				size_t Index_CF32off, Index_SF32off, Index_CI64off, Index_SI64off;
				
			};
			CompartmentGroupingImplementation comp_group_impl;
			
			// for Cable Equation solvers
			struct CableSolverDefinition{
				
				SimulatorConfig::CableEquationSolver type;
				
				// default structure from NeuroML, till unique cells LATER
				std::vector< Int > BwdEuler_OrderList;
				std::vector< Int > BwdEuler_ParentList;
				std::vector< Real > BwdEuler_InvRCDiagonal;
				
			}cable_solver;
			
			
			
			
			struct CableSolverImplementation{
				// for tables, for each compartment
				size_t Index_BwdEuler_OrderList; // table
				size_t Index_BwdEuler_ParentList; // table
				size_t Index_BwdEuler_InvRCDiagonal; // table
				
				size_t Index_BwdEuler_WorkDiagonal; // table
			};
			CableSolverImplementation cable_solver_implementation;
			
			
			struct IonChannelDistributionInstance{
				
				Int ion_species;
				Int ion_channel;
				ChannelDistribution::Type type;
				
				Real conductivity;
				Real erev; // for Fixed and Population
				Real vshift; //for vshift
				Real permeability; // for GHK1
				Int number; // for Population
				
				IonChannelDistributionInstance(){
					
				}
			};
			struct IonSpeciesDistributionInstance{
				
				//Int ion_species;
				Int conc_model_seq;
				
				Real initialConcentration; // used with core types
				Real initialExtConcentration;
				
				IonSpeciesDistributionInstance(){
					
				}
			};
			// TODO place somewhere more appropriate
			// except for data values, it's the same between structurally-identical compartments
			// so they can be simulated with the exact same code
			struct CompartmentDefinition{
				// could convert into proxy object, to allow Struct-of-Arrays transformation, for both memory and time efficiency, LATER
				Real V0, Vt, AxialResistance, Capacitance;
				
				
				std::vector<IonChannelDistributionInstance> ionchans;
				std::map< Int, IonSpeciesDistributionInstance > ions;
				
				std::vector<int32_t> adjacent_compartments;
				
				IdListRle input_types;
				IdListRle synaptic_component_types;
				
				bool spike_output;
				
				// TODO hash
				// TODO perhaps create a structure subset, that is enough to generate the code and differentiate the codes
				
			};
			std::vector< CompartmentDefinition > comp_definitions;
			
			// indices to retrieve which table is which
			// mirror the structure of the cell's components in the NeuroML definition
			
			struct IonChannelDistImplementation{
				struct SubGate{
					Int Index_Q10; // when used
					Int Index_Q10_BaseTemp; // when used
					
					Int Index_Q; // if KS gate transition
					
					ComponentSubSignature alpha_component; // which component instance?
					ComponentSubSignature  beta_component; // which component instance? 
					ComponentSubSignature   tau_component; // which component instance? 
					ComponentSubSignature   inf_component; // which component instance? also used for instantaneous gates
					SubGate(){
						Index_Q = -1;
					}
				};
				struct PerGate{
					Int Index_Q10; // when used
					Int Index_Q10_BaseTemp; // when used
					
					Int Index_Q; // -1 if composite, KS gate etc.
					
					ComponentSubSignature alpha_component; // which component instance?
					ComponentSubSignature  beta_component; // which component instance? 
					ComponentSubSignature   tau_component; // which component instance? 
					ComponentSubSignature   inf_component; // which component instance? also used for instantaneous gates, and lems gates too!
					
					// duplicated from SubGate, in case an important difference appears
					
					std::vector<SubGate> subgates;
					std::vector<SubGate> transitions;
					
					PerGate(){
						Index_Q = -1;
					}
				};
				
				struct ConductanceScaling{
					Int Index_Q10; // when used
					Int Index_Q10_BaseTemp; // when used
					
					ComponentSubSignature scaling_component;
				};
				
				// if native type
					ConductanceScaling conductance_scaling;
					std::vector<PerGate> per_gate;
				// if LEMS component instead
					ComponentSubSignature channel_component; 
			};
			
			struct IonSpeciesDistImplementation{
				// type is defined at the per-compartment distribution array, internally
				
				// for every distribution instance
				size_t Index_InitIntra; // constant
				size_t Index_InitExtra; // constant
				
				// when it's a native type
				size_t Index_Intra; // state
				size_t Index_Extra; // state
				size_t Index_RestConc; // constant
				size_t Index_DecayTau; // constant
				size_t Index_Shellthickness_Or_RhoFactor; // constant
				
				// when it's a component
				ComponentSubSignature component; 
				
				IonSpeciesDistImplementation(){
					Index_Intra = -1;
					Index_Extra = -1;
				}
			};
			
			struct CompartmentImplementation{
				std::map< Int, IonSpeciesDistImplementation > concentration; // per ion species
				std::vector<IonChannelDistImplementation> channel; // per distribution
				std::map< Int, InputImplementation > input; // per input type id_id
				std::map< Int, SynapticComponentImplementation > synapse; // per synapse type id_id
				
				SpikeSendingImplementation spiker; // just one for the present model of spike-sending
				
				ptrdiff_t Index_AdjComp; // TODO turn into full implementation, perhaps Table_?
				
				CompartmentImplementation(){
					Index_AdjComp = -1;
				}
			};
			std::vector<CompartmentImplementation> comp_implementations; // per segment
			
			
			size_t GetVoltageStatevarIndex( Int comp_seq ) const {
				// alternate memory layouts LATER
				return Index_Voltages + comp_seq;
			}
		};
		// just for convenience
		typedef PhysicalCell::IonChannelDistributionInstance IonChannelDistributionInstance;
		typedef PhysicalCell::IonSpeciesDistributionInstance IonSpeciesDistributionInstance;
		typedef PhysicalCell::IonChannelDistImplementation IonChannelDistImplementation;
		typedef PhysicalCell::IonSpeciesDistImplementation IonSpeciesDistImplementation;
		typedef PhysicalCell::CompartmentGrouping CompartmentGrouping;
		typedef PhysicalCell::CompartmentImplementation CompartmentImplementation;
		typedef PhysicalCell::CompartmentDefinition CompartmentDefinition;
		
		PhysicalCell physical_cell;
		
		// The mapping of properties to offsets, for artificial cells (aka point neurons)
		struct ArtificialCell{
			// <attachment>s common in networks, regardless of cell type(type)
			std::map< Int, InputImplementation > input; // per input type id_id
			std::map< Int, SynapticComponentImplementation > synapse; // per synapse type id_id
			
			SpikeSendingImplementation spiker; // just one for the present model of spike-sending
			
			// if a LEMS component
			//ComponentValueSignature input_component_values;
			ComponentSubSignature component;
			
			ptrdiff_t Index_Statevar_Voltage; // if present
			
			// just so, handle it better LATER
			InputImplementation inpimpl;
			
			ArtificialCell(){
				Index_Statevar_Voltage = -1;
			}
		};
		ArtificialCell artificial_cell;
		
		// the mapping of what is common between cell implementations
		// will refactor better LATER, as converting individual cells to blocks goes
		struct CommonInCell{
			
			RngImplementation cell_rng_seed;
		};
		CommonInCell common_in_cell;
		
		
		// the actual signature of the component
		//TODO class/instance constants
		struct WorkItemDataSignature{
			
			// XXX find out why AppendToVector generic path has such spectacular slowdown,
			// 	Table_F32 vs. std::vector<float>. Could it be the allocator's fault, STL non-inlining's fault, or what else?
			// std::vector<float> state, constants;
			RawTables::Table_F32 state, constants;
			
			std::vector<TableInfo> tables_const_f32, tables_const_i64, tables_state_f32, tables_state_i64;
			std::unordered_map<size_t, std::string> constants_names;
			std::unordered_map<size_t, std::string> state_names;
			
			// TODO do something with these, default constants for tables??
			std::vector<Real> prototype_const;
			std::vector<Real> prototype_state;
			// LATER reverse indices, somehow
			
			// XXX call counter should be replaced by structured indices for each expression in use
			// to avoid having the symbolically same invocation produce different results  !
			// on the other hand, this bug may occur when:
			// 	a component is exposed twice
			// 	and its input is important (indirectly exposed) in both cases.
			// Right now this only happens with concentration models, which typically don't use RNG for exposures,
			// 	so this is not a fix-me outright
			Int random_call_counter; 
			
			void Append( const WorkItemDataSignature &rhs ){
				
				ptrdiff_t conoff = constants.size();
				for( auto keyval : rhs.constants_names ){
					constants_names[ keyval.first + conoff ] = keyval.second;
				}
				ptrdiff_t staoff = state.size();
				for( auto keyval : rhs.state_names ){
					state_names[ keyval.first + staoff ] = keyval.second;
				}
				
				AppendToVector( constants        , rhs.constants        );
				AppendToVector( state            , rhs.state            );
				AppendToVector( tables_const_f32 , rhs.tables_const_f32 );
				AppendToVector( tables_state_f32 , rhs.tables_state_f32 );
				AppendToVector( tables_const_i64 , rhs.tables_const_i64 );
				AppendToVector( tables_state_i64 , rhs.tables_state_i64 );
				
				AppendToVector( prototype_const  , rhs.prototype_const  );
				AppendToVector( prototype_state  , rhs.prototype_state  );
				
				random_call_counter += rhs.random_call_counter;
			}
			
			WorkItemDataSignature(){
				random_call_counter = 0;
			}
		};
		
		WorkItemDataSignature cell_wig;
		
		std::string code;
		IterationCallback callback;
		std::string name;
	};
	std::vector<CellInternalSignature> cell_sigs;
	// Allocate empty signatures for neuron types, and start filling them in in multiple pases
	cell_sigs.resize(cell_types.contents.size());
	
	//------------------>  Determine just how neurons are discretised. This will aid in determining which synapses are located on which compartments.
	printf("Analyzing cell morphologies...\n");
	for(size_t cell_seq = 0; cell_seq < cell_types.contents.size(); cell_seq++){
		const auto &cell_type = cell_types.contents[cell_seq];
		auto &sig = cell_sigs[cell_seq];
		
		if( cell_type.type == CellType::PHYSICAL ){
			const PhysicalCell &cell = cell_type.physical;
			const Morphology &morph = morphologies.get(cell.morphology);
			
			const bool debug_log_discretisation = (config.debug && 0);
			
			// Before determining the values of compartments, it would be convenient to lay out their symbolic layout on an array.
			// Fortunately, we can determine this for the NeuroML 'cables' description, before actually chopping the neurite segments;
			// though this might not be generalisable for other algorithms.
			// The invariant of the array representation of compartments is the same as the array representation of segments: that parents must come before children.
			// Therefore the cables may appear in a different order than the one they were specified in the original NeuroML description, since segment groups may appear in mostly any order.
			
			// Determine which segments belong to which 'cable' group.
			std::vector<Int> cable_per_segment(morph.segments.size(), -1); //group_seq of cable for each segment
			
			const bool config_use_cable_discretization = 1;
			if( config_use_cable_discretization ){
			for( Int group_seq = 0; group_seq < (Int) morph.segment_groups.size(); group_seq++){
				const auto &group = morph.segment_groups[group_seq];
				
				if( !group.is_cable ) continue;
				
				std::vector<Int> segs_in_cable = group.list.toArray();
				for( Int seg_seq : segs_in_cable ){
					// perhaps the segment is already included in another group?
					Int previous_group_seq = cable_per_segment[seg_seq];
					if(!(previous_group_seq < 0)){
						// cables should not overlap!!
						fprintf(stderr, "internal error: segment groups %ld and %ld should be non-overlapping, yet they overlap on segment %ld", group_seq, previous_group_seq, morph.lookupNmlId(seg_seq));
						return false;
					}
					cable_per_segment[seg_seq] = group_seq;
				}
				
			}
			}
			
			// Assign compartments to the parts of the neuron. Each cable is assigned a continuous span of compartments, each non-cable segment is assigned its own compartment.
			std::vector<Int> free_segment_to_compartment(morph.segments.size(), -1);
			std::vector<Int> cable_to_compartment_offset(morph.segment_groups.size(), -1);
			// TODO clarify: nseg is the number of *compartments* per cable
			std::vector<Int> nseg_per_cable(morph.segment_groups.size(), -1);
			
			// LATER d_lambda rule perhaps?
			Int number_of_compartments = 0;
			for( Int seg_seq = 0; seg_seq < (Int)morph.segments.size(); seg_seq++ ){
				Int cable_seq = cable_per_segment[seg_seq];
				
				if( cable_seq >= 0 ){
					// a cable
					const Morphology::SegmentGroup &group = morph.segment_groups[cable_seq];
					
					assert(group.is_cable);
					// allocate only if not the cable has not been allocated already, by a previous segment
					if(!(cable_to_compartment_offset[cable_seq] < 0)) continue;
					
					Int cable_nseg = group.cable_nseg;
					if( cable_nseg < 0 ) cable_nseg = 1; // by default, TODO add a shorthand getter in model
					nseg_per_cable[cable_seq] = cable_nseg;
					
					cable_to_compartment_offset[cable_seq] = number_of_compartments;
					number_of_compartments += cable_nseg;
				}
				else{
					// a free-standing segment
					// since NML segments are being sweeped, this free-standing segment should be encountered only once in this loop
					assert(free_segment_to_compartment[seg_seq] < 0);
					free_segment_to_compartment[seg_seq] = number_of_compartments;
					number_of_compartments++;
				}
			}
			if(debug_log_discretisation){
				printf("Number of compartments in cell: %ld\n", number_of_compartments);
				printf("Compartments of Free-standing segments:\n");
				for( Int seg_seq = 0; seg_seq < (Int)morph.segments.size(); seg_seq++ ){
					if( free_segment_to_compartment[seg_seq] >= 0 ){
						printf("%ld ",seg_seq);
					}
				}
				printf("\n");
				printf("Compartments of Cables:\n");
				for( Int group_seq = 0; group_seq < (Int)morph.segment_groups.size(); group_seq++ ){
					Int first_comp = cable_to_compartment_offset[group_seq];
					if( first_comp < 0 ) continue;
					
					printf("%ld - %ld ", first_comp, nseg_per_cable[group_seq] - 1 );
				}
				printf("\n");
			}
			auto &pig = sig.physical_cell;
			auto &comp_disc  = pig.compartment_discretization;
			auto &segment_to_compartment_seq                       = comp_disc.segment_to_compartment_seq                      ;
			auto &segment_to_compartment_fractionAlong             = comp_disc.segment_to_compartment_fractionAlong            ;
			auto &compartment_to_segment_seq                       = comp_disc.compartment_to_segment_seq                      ;
			auto &compartment_to_first_segment_start_fractionAlong = comp_disc.compartment_to_first_segment_start_fractionAlong;
			auto &compartment_to_last_segment_end_fractionAlong    = comp_disc.compartment_to_last_segment_end_fractionAlong   ;
			
			
			comp_disc.number_of_compartments = number_of_compartments;
			// Prepare the maps from compartments to cells, by allocating them to the right size
			// for each compartment, keep a list of the segments it contains, and fractionAlong for beginning and end.
			// (Implicitly, each segment in between is included until its end and from its start.)
			segment_to_compartment_seq.resize( morph.segments.size() ); 
			segment_to_compartment_fractionAlong.resize( morph.segments.size() );
			
			compartment_to_segment_seq.resize( number_of_compartments ); 
			compartment_to_first_segment_start_fractionAlong.assign( number_of_compartments, -1 ); 
			compartment_to_last_segment_end_fractionAlong.assign( number_of_compartments, -1 ); 
			
			// Process the structure of cables: slice them into compartments of equal length.
			// later on, the density mechanism specs have to be applied on segments anyway,
			// 	so the discretisation can be pre-scribed without evaluating any physical comp properties (other than path length) in advance.
			for(Int group_seq = 0; group_seq < (Int)morph.segment_groups.size(); group_seq++){
				const Morphology::SegmentGroup &group = morph.segment_groups[group_seq];
				
				if(!(group.is_cable && config_use_cable_discretization)) continue;
				
				auto cable_nseg = nseg_per_cable[group_seq];
				
				// Iterate the segments, chop them into the prescribed number of compartments.
				// First, count the length of the cable
				std::vector<Real> segment_length_over_cable(group.list.Count(),NAN);
				Real running_total_segment_length = 0;
				
				// the sorted list of seg_seq 
				std::vector<Int> segs_in_cable = group.list.toArray(); 
				for(Int i = 0; i < (Int)segs_in_cable.size(); i++){
					const auto &seg = morph.segments.atSeq(segs_in_cable[i]); 
					if( i > 0 ){
						if(seg.parent != segs_in_cable[i-1]){
							fprintf(stderr, "internal error: group is not a proper unbranched cable: segment %ld 's parent is %ld when it should be %ld",
								morph.segments.getId(segs_in_cable[i]), morph.segments.getId(seg.parent), morph.segments.getId(segs_in_cable[i-1]));
							return false;
						}
						// check that fractionAlong = 1 (attached to end of parent) for a properly unbranched cable
						if(seg.parent >= 0 && seg.fractionAlong != 1){
							fprintf(stderr, "internal error: group is not a proper unbranched cable: segment %ld is not attached to end of parent but on fractionAlong = %.17g instead",
								morph.segments.getId(segs_in_cable[i]), seg.fractionAlong);
							return false;
						}
					}
					running_total_segment_length += DistancePt3d( seg.proximal, seg.distal );
					segment_length_over_cable[i] = running_total_segment_length;
					
				}
				// running_total_segment_length has total length of the cable by now
				const auto total_cable_length = running_total_segment_length;
				
				if(debug_log_discretisation){
					printf("segcumlen\n");
					for(int i = 0; i < (Int)segs_in_cable.size(); i++){
						printf("%f\n", segment_length_over_cable[i]);
					}
				}
				
				// now slide along the cable, and chop/merge segments into the new compartments
				Int segment_under_compartment = 0; // invariant: after each iteration, segment_under_compartment points to the segment over the farthest edge of the last processed compartment
				// therefore end_compartment_path_length <= segment_length_over_cable[segment_under_compartment]
				
				// NB this variable is not used since compartment borders are distributed evenly along path length of the cable. 
				// It might be useful in case the distribution is not even LATER
				// std::vector<Real> compartment_length_over_cable(cable_nseg,NAN);
				// TODO special case where soma is spherical and yet participating in a cable! Keep the soma spherical and distribute nseg - 1?
				// current implementation will simply not account for the segment as it has zero length
				for(Int comp_in_cable = 0; comp_in_cable < cable_nseg; comp_in_cable++){
					
					// The borders of path_length that this compartment covers
					// The mapping is fixed and linear over path_length, for now
					Real start_compartment_path_length = ((comp_in_cable+0)/Real(cable_nseg)) * total_cable_length; // Divide by nseg before multiplying other factors, to get a perfect 1.0 for the last seg
					Real end_compartment_path_length = ((comp_in_cable+1)/Real(cable_nseg)) * total_cable_length;
					
					Int comp_seq = cable_to_compartment_offset[group_seq] + comp_in_cable;
					
					if(debug_log_discretisation) printf("cable group_seq %d comp %d, comp_seq %d begin\n", (int) group_seq, (int)comp_in_cable, (int) comp_seq);
					
					
					// Gather the slices of NML segments
					// Since all said properties are directly summable from their component parts, this is an inline method to avoid extra complexity, for now.
					auto AddSegmentSliceToCompartment = [ 
						&segment_to_compartment_seq, &segment_to_compartment_fractionAlong,
						&compartment_to_segment_seq, &compartment_to_first_segment_start_fractionAlong, &compartment_to_last_segment_end_fractionAlong
					]( 
						Int comp_seq, Int seg_seq, Real segment_fractionAlong_start, Real segment_fractionAlong_end 
					){
						assert(0 <= segment_fractionAlong_start && segment_fractionAlong_start <= segment_fractionAlong_end && segment_fractionAlong_end <= 1);
						
						// for each segment, keep a list of the compartments it spans, and up to which fractionAlong it spans them.
						// for each compartment, keep a list of the segments it contains, and fractionAlong for beginning and end.
						
						// if it is the first segment being added to the compartment, set the fractionAlong
						if( compartment_to_segment_seq[comp_seq].empty() ){
							compartment_to_first_segment_start_fractionAlong[comp_seq] = segment_fractionAlong_start;
						}
						// set the fractionAlong of last segment incrementally, to be updated by following segments being attached.
						compartment_to_last_segment_end_fractionAlong[comp_seq] = segment_fractionAlong_end;
						// maybe setting the first and last fractionAlong can be done explicitly in the per-cable chopping loop, but this is less fragile.
						
						
						compartment_to_segment_seq[comp_seq].push_back(seg_seq);
						
						segment_to_compartment_seq[seg_seq].push_back(comp_seq);
						segment_to_compartment_fractionAlong[seg_seq].push_back(segment_fractionAlong_end);
						
					};
					
					// helpers for converting path length to fractionAlong, whet else is there to say, with a fallback for zero length segments
					// TODO make it more elegant
					auto PathLengthToFractionAlong_OrZero = [](Real start, Real end, Real query){
						Real fractionAlong = (query - start)/(end - start);
						if( abs(end - start) < 0.001 ){
							// in case the diameter changes in an abruptly short span of path length, fix segment_fractionAlong to a valid and quite accurate value
							// NB: this remains ambiguous! fractionAlong could be 1 as well, therefore this must be checked beforehand or this routine should be called on a codepath where 1 is explicitly used
							fractionAlong = 0;
						}
						return fractionAlong;
					};
					auto SegmentQueryToFractionAlong_OrZero = [PathLengthToFractionAlong_OrZero, &segment_length_over_cable](Int seg_idx, Real query_path_along){
						Real start_segment_path_length = 0;
						if( seg_idx > 0 ) start_segment_path_length = segment_length_over_cable[seg_idx-1];
						Real end_segment_path_length = segment_length_over_cable[seg_idx];
						
						return PathLengthToFractionAlong_OrZero(start_segment_path_length, end_segment_path_length, query_path_along);
					};
					
					// The path_length that has been accumulated in the compartment so far
					Real index_compartment_path_length = start_compartment_path_length;
					
					// merge all segments that end up to end_compartment_path_length
					while( 
						segment_under_compartment < (Int)segment_length_over_cable.size() 
						&& segment_length_over_cable[segment_under_compartment] <= end_compartment_path_length 
					){
						// Integrate frustum from index_compartment_path_length to end of current segment.
						
						Real segment_fractionAlong = SegmentQueryToFractionAlong_OrZero(segment_under_compartment, index_compartment_path_length); // NB: returns zero for empty segment
						if( abs(segment_fractionAlong - 1) < 0.0001 ){
							// no need to integrate this minuscule part, skip
							//would add abs(index_compartment_path_length - segment_length_over_cable[segment_under_compartment]) < 0.001 condition but have to account for zero length segments
						}
						else{
							auto seg_seq = segs_in_cable[segment_under_compartment];
							if(debug_log_discretisation) printf("complete seg %ld (%d) fractionAlong %f->end\n", segment_under_compartment, (int) seg_seq, segment_fractionAlong);
							AddSegmentSliceToCompartment(comp_seq, seg_seq, segment_fractionAlong, 1);
						}
						
						index_compartment_path_length = segment_length_over_cable[segment_under_compartment];
						segment_under_compartment++;
					}
					// If the compartments spans a part of the, 
					if( segment_under_compartment < (Int)segment_length_over_cable.size() 
						&& end_compartment_path_length < segment_length_over_cable[segment_under_compartment] 
					){
						// Integrate frustum from index_compartment_path_length to end of current compartment, on current segment.
						Real segment_fractionAlong_start = SegmentQueryToFractionAlong_OrZero(segment_under_compartment, index_compartment_path_length);
						Real segment_fractionAlong_end = SegmentQueryToFractionAlong_OrZero(segment_under_compartment, end_compartment_path_length);
						// add an exception to include the whole segment to the end, if it is the last compartment of the cable. To avoid letting this end be unaccounted for.
						if( comp_in_cable + 1 == cable_nseg ){
							// segment_fractionAlong_end = 1;
						}
						
						
						if( abs(segment_fractionAlong_end - segment_fractionAlong_start) < 0.0001 ){
							// no need to integrate this minuscule part, skip
						}
						else{
							if(debug_log_discretisation) printf("partial seg %ld fractionAlong %f->%f\n", segment_under_compartment, segment_fractionAlong_start, segment_fractionAlong_end);
							AddSegmentSliceToCompartment(comp_seq, segs_in_cable[segment_under_compartment], segment_fractionAlong_start, segment_fractionAlong_end);
						}
						
						// compartment_length_over_cable[comp_in_cable] = end_compartment_path_length;
					}
				}
				
			}
			
			// Define the segments not in cables as individual compartments.
			// Perhaps save memory somehow, in this case LATER
			for( Int seg_seq = 0; seg_seq < (Int)morph.segments.size(); seg_seq++ ){
				Int comp_seq = free_segment_to_compartment[seg_seq];
				if(!( comp_seq >= 0 )) continue;
				
				segment_to_compartment_seq[seg_seq].push_back(comp_seq);
				segment_to_compartment_fractionAlong[seg_seq].push_back(1);
				
				compartment_to_segment_seq[comp_seq].push_back(seg_seq);
				compartment_to_first_segment_start_fractionAlong[comp_seq] = 0;
				compartment_to_last_segment_end_fractionAlong[comp_seq] = 1;
				
			}
			
			// Resolving the (parent-child) connections between segments to the corresponding compartments
			// is done here, it is more convenient since the cables are processed in this pass over cell types.
			// This could be rolled in the preceding loops (per cable and per segment) but it would complicate that code further.
			// (Also by requiring iterating over cables to follow a seg_seq_compatible order so that prent cable is already resolved)
			// Another alternative is to ierate through each segment, linking compartments in series along each,
			// 	and linking compartments that the segment starts with since fractionAlong=0, to their parents. But why make this more difficult than it has to be
			// {
			auto &tree_parent_per_compartment = comp_disc.tree_parent_per_compartment;
			tree_parent_per_compartment.assign( number_of_compartments, -1 );
			// std::vector<bool> cable_processed( morph.segment_groups.size(), false );
			
			// Resolve parents for free segments.
			for( Int seg_seq = 0; seg_seq < (Int)morph.segments.size(); seg_seq++ ){
				Int comp_seq = free_segment_to_compartment[seg_seq];
				if( comp_seq < 0 ) continue;
				
				//the segment maps to exactly one compartment, map it
				const auto &seg = morph.segments.atSeq(seg_seq); 
				if( seg.parent >= 0 ){
					tree_parent_per_compartment[comp_seq] = comp_disc.GetCompartmentForSegmentLocation( seg.parent, seg.fractionAlong );
				}
				else tree_parent_per_compartment[comp_seq] = -1;
				
			}
			// Resolve parents for cables.
			for(Int group_seq = 0; group_seq < (Int)morph.segment_groups.size(); group_seq++){
				const Morphology::SegmentGroup &group = morph.segment_groups[group_seq];
				
				if(!(group.is_cable && config_use_cable_discretization)) continue;
				
				
				auto cable_nseg = nseg_per_cable[group_seq];
				assert( cable_nseg >= 1 );
				
				auto first_comp_on_cable = cable_to_compartment_offset[group_seq];
				
				// For all compartments succeeding the first one, their parent property is obvious: i is the previous comp_seq in order obviously to their predecessors
				for( int comp_seq = first_comp_on_cable + 1; comp_seq < first_comp_on_cable + cable_nseg; comp_seq++){
					tree_parent_per_compartment[comp_seq] = comp_seq - 1;
				}
				// For the first compartment, use the parent of the first seg of the fisrt compartment
				// (the first seg should, of course, be included in the first compartment starting from fractionAlong = 0)
				Int first_seg_on_cable = compartment_to_segment_seq[first_comp_on_cable][0];
				
				assert(compartment_to_first_segment_start_fractionAlong[first_comp_on_cable] == 0);
				
				const auto &seg = morph.segments.atSeq(first_seg_on_cable); 
				if( seg.parent >= 0 ){
					tree_parent_per_compartment[first_comp_on_cable] = comp_disc.GetCompartmentForSegmentLocation( seg.parent, seg.fractionAlong );
				}
				else tree_parent_per_compartment[first_comp_on_cable] = -1;
				
			}
			
			if(debug_log_discretisation){
				printf("Mapping:\n");
				for( Int seg_seq = 0; seg_seq < (Int)morph.segments.size(); seg_seq++ ){
					printf("Seg %d:", (int) seg_seq);
					for( Int j = 0; j < (Int) segment_to_compartment_seq[seg_seq].size(); j++ ){
						printf(" %d (%g)", (int) segment_to_compartment_seq[seg_seq][j], segment_to_compartment_fractionAlong[seg_seq][j]);
					}
					printf("\n");
					
				}
				for( Int comp_seq = 0; comp_seq < (Int)compartment_to_segment_seq.size(); comp_seq++ ){
					printf("Comp %d:", (int) comp_seq);
					printf(" (%g)", compartment_to_first_segment_start_fractionAlong[comp_seq]);
					for( Int j = 0; j < (int) compartment_to_segment_seq[comp_seq].size(); j++ ){
						printf(" %d", (int) compartment_to_segment_seq[comp_seq][j]);
					}
					printf(" (%g)", compartment_to_last_segment_end_fractionAlong[comp_seq]);
					printf("\n");
				}
				
				printf("Parents:\n");
				for( Int comp_seq = 0; comp_seq < number_of_compartments; comp_seq++ ) 	printf("\t%d", (int) comp_seq);	
				printf("\n");
				for( Int comp_seq = 0; comp_seq < number_of_compartments; comp_seq++ ) 	printf("\t%d", (int) tree_parent_per_compartment[comp_seq]);	
				printf("\n");
			}
			
		}
		else if( cell_type.type == CellType::ARTIFICIAL ){
			// no segments no worries
		}
		else{
			assert(false);
		}
	}
	// assert(false);
	printf("Analyzing connectivity...\n");
	//------------------>  Scan inputs, to aid cell type analysis
	// Current approach: For each cell type, aggregate a set of possible incident synapse or input mechanisms.
	// This is far from the only appropriate analysis method. 
	// 	For example, it could be smarter to discriminate between populations of the same cell types, or subsets of individual instances.
	// 	Or even forgoing grouping by cell, instead grouping compartments with the same sets of attachment types.
	// The current approach is targeting the existing code generation approach of making one code block for each cell type.
	// per cell, per segment
	auto GetInputIdId = [ &input_sources ]( Int input_seq ){
		const InputSource &inp = input_sources.get( input_seq );
		Int id_id;
		if(inp.type == InputSource::Type::COMPONENT) id_id = input_seq; // syn.component; TODO
		if(
			inp.type == InputSource::Type::TIMED_SYNAPTIC
			|| inp.type == InputSource::Type::POISSON_SYNAPSE
			|| inp.type == InputSource::Type::POISSON_SYNAPSE_TRANSIENT
		) id_id = input_seq; // structural differences
		else if( inp.component.ok() ){
			id_id = input_seq; // implement the devious trick here XXX
		}
		else id_id = Int(inp.type) - InputSource::Type::MAX;
		return id_id;
	};
	
	// NOTE this is just a helper for this sort of mapping to compartments (includes mapping all to compartment zero for artificial cells).
	// TODO either roll it into CellInternalSignature or restrict its scope.
	auto GetCompSeqForCell = []( const CellType &cell_type, const CellInternalSignature &sig,  Int seg_seq, Real fractionAlong ){
		if( cell_type.type == CellType::Type::PHYSICAL ){
			return sig.physical_cell.compartment_discretization.GetCompartmentForSegmentLocation(seg_seq, fractionAlong);
		}
		else if( cell_type.type == CellType::Type::ARTIFICIAL ){
			return 0;
		}
		else{
			assert(false);
			return -1;
		}
	};
	
	std::vector< std::map<Int, IdListRle> > input_types_per_cell_( cell_types.contents.size() );
	std::vector< std::map<Int, IdListRle> > input_types_per_cell_per_compartment( cell_types.contents.size() );
	for(size_t inp_seq = 0; inp_seq < net.inputs.size(); inp_seq++){
		const auto &inp = net.inputs[inp_seq];
		
		const auto &pop = net.populations.get(inp.population);
		
		const Int comp_seq = GetCompSeqForCell( cell_types.get(pop.component_cell), cell_sigs[pop.component_cell], inp.segment, inp.fractionAlong );
		
		// LATER could specialize code to elide weight
		int id_id = GetInputIdId( inp.component_type );
		// TODO perhaps discriminate between same-type input types, to re-use globals LATER
		input_types_per_cell_per_compartment[pop.component_cell][comp_seq].Addd(id_id); //TODO when work unit is compartment, probably after morpho analysis
		
	}
	// and normalize input_types_per_cell_per_compartment lists after that
	for(size_t cell_seq = 0; cell_seq < cell_types.contents.size(); cell_seq++){
		for(auto keyval : input_types_per_cell_per_compartment[cell_seq]){
			keyval.second.Compact();
		}
	}
	
	//------------------>  Scan synaptic projections, to aid cell type analysis
	// TODO what about per-compartment splitting? //perhaps generate comparmtent-internal constants and then combine with synapses to generate code
	
	// Gather all pre- and post- synaptic components applicable on a cell type
	// todo name bijection
	// per cell, set of segments
	std::vector< std::set<Int> > spiking_outputs_per_cell_( cell_types.contents.size() );
	std::vector< std::set<Int> > spiking_outputs_per_cell_per_compartment( cell_types.contents.size() );
	
	// per cell, per segment, list of components
	// for aggregation of tabular types TODO refactor
	// clashes with deduplication of constants
	auto GetSynapseIdId = [ &synaptic_components ]( Int syncomp_seq ){
		const SynapticComponent &syn = synaptic_components.get( syncomp_seq );
		Int id_id;
		if(syn.type == SynapticComponent::Type::COMPONENT) id_id = syncomp_seq; // syn.component; TODO
		// special cases requiring special treatment, despite using LEMS component
		else if( syn.type == SynapticComponent::BLOCKING_PLASTIC ) id_id = syncomp_seq;
		else if( syn.component.ok() ){
			id_id = syncomp_seq; // implement the devious trick here XXX
		}
		else id_id = Int(syn.type) - SynapticComponent::Type::MAX;
		return id_id;
	};
	std::vector< std::map<Int, IdListRle> > synaptic_component_types_per_cell_( cell_types.contents.size() );
	std::vector< std::map<Int, IdListRle> > synaptic_component_types_per_cell_per_compartment( cell_types.contents.size() );
	
	for(size_t proj_seq = 0; proj_seq < net.projections.contents.size(); proj_seq++){
		const auto &proj = net.projections.get(proj_seq);
		const auto &prepop = net.populations.get(proj.presynapticPopulation);
		const auto &postpop = net.populations.get(proj.postsynapticPopulation);
		
		
		for(const auto &conn : proj.connections.contents){
			
			const Int pre_comp_seq  = GetCompSeqForCell( cell_types.get( prepop.component_cell  ), cell_sigs[prepop.component_cell ], conn.preSegment , conn.preFractionAlong  );
			const Int post_comp_seq = GetCompSeqForCell( cell_types.get( postpop.component_cell ), cell_sigs[postpop.component_cell], conn.postSegment, conn.postFractionAlong );
		
			
			// If a syn.component needs Vpeer, make indices for Vpeer
			// If a syn.component needs spike, make indices for peer to send spikes
			
			// LATER could specialize code to elide weight
			// LATER could specialize code to elide delay
			
			if(conn.type == Network::Projection::Connection::SPIKING){
				// TODO validate conditions:
				// Vt must be specified for spike producing compartments (also in NeuroML API)
				// pre-synaptic must have a spike output port (in LEMS)
				// post-synaptic must have a spike input (in LEMS)
				
				int id_id = GetSynapseIdId( conn.synapse );
				
				spiking_outputs_per_cell_per_compartment[prepop.component_cell ].insert(pre_comp_seq );
				synaptic_component_types_per_cell_per_compartment[postpop.component_cell][post_comp_seq].Addd(id_id);
			}
			else if(conn.type == Network::Projection::Connection::ELECTRICAL){
				// same behaviour for pre-and post-synaptic
				// Vpeer is required, and it always exists for physical compartments
				
				//printf("yeo %ld %ld %ld\n\n",conn.synapse, prepop.component_cell, pre_comp_seq);
				int id_id = GetSynapseIdId( conn.synapse );
				
				synaptic_component_types_per_cell_per_compartment[prepop.component_cell ][pre_comp_seq ].Addd(id_id);
				synaptic_component_types_per_cell_per_compartment[postpop.component_cell][post_comp_seq].Addd(id_id);
			}
			else if(conn.type == Network::Projection::Connection::CONTINUOUS){
				
				// could use either, TODO refactor attachment to reflect this
				// this is the most general case, perhaps keep it? TODO
				int id_id_pre  = GetSynapseIdId( conn.continuous.preComponent  );
				int id_id_post = GetSynapseIdId( conn.continuous.postComponent );
				
				if( synaptic_components.get( conn.continuous.postComponent ).HasSpikeIn( component_types ) ){
					spiking_outputs_per_cell_per_compartment[prepop.component_cell ].insert(pre_comp_seq );
				}
				
				if( synaptic_components.get( conn.continuous.preComponent ).HasSpikeIn( component_types ) ){
					spiking_outputs_per_cell_per_compartment[postpop.component_cell].insert(post_comp_seq);
				}
				
				synaptic_component_types_per_cell_per_compartment[prepop.component_cell ][pre_comp_seq ].Addd(id_id_pre );
				synaptic_component_types_per_cell_per_compartment[postpop.component_cell][post_comp_seq].Addd(id_id_post);
			}
			else{
				// TODO internal error
				return false;
			}
		}
	}
	// and normalize synaptic_component_types_per_cell_per_compartment lists after that
	for(size_t cell_seq = 0; cell_seq < cell_types.contents.size(); cell_seq++){
		for(auto keyval : synaptic_component_types_per_cell_per_compartment[cell_seq]){
			keyval.second.Compact();
		}
	}
	// TODO also scan event logging for emitters
	// LATER perhaps also scan trajectory loggers too, for MPI or what? just preserve the dependencies
	
	//------------------>  Analyze cell types
	
	// also generate the unit conversion suffixes to keep engine units consistent
	// a consistent selection of engine units should need none, but let it be configurable to permit technicalities (and automatically prevent bugs)
	struct Convert{
		static std::string Suffix(const ScaleEntry &scale){
			std::string ret;
			char tmps[50];
			if(scale.scale != 1){
				sprintf(tmps, " * %.17g", scale.scale); // maintain full precision in decimal output; exactly the same value will be passed to compiler
				ret += tmps;
			}
			if(scale.pow_of_10 != 0){
				sprintf(tmps, " * 1e%df", scale.pow_of_10);
				ret += tmps;
			}
			if(scale.offset != 0){
				sprintf(tmps, " + %.17g", scale.offset);
				ret += tmps;
			}
			return ret;
		}
	};
	
	// append stuff to work item, whether one-off or multi-instance
	struct ISignatureAppender{
		
		virtual size_t Constant( Real default_value, const std::string &for_what ) const = 0;
		virtual size_t StateVariable( Real default_value, const std::string &for_what ) const = 0;
		
		virtual std::string ReferTo_Const( size_t index ) const = 0;
		virtual std::string ReferTo_State( size_t index ) const = 0;
		virtual std::string ReferTo_StateNext( size_t index ) const = 0;
	};
	/*
	LEMS component implementation:
		- Declare the constants and state variables
		- Calculate derived variables whenever required
		- Update dynamics - requires definitions
		- Assign initial state - requires definitions
		- Invoke calculations from engine itself ? - requires symbolic value closure
		- logging of derived vars - invoke from engine or expose as state variable
		
		LEMS component with in port "IN" and out port "OUT" :
		
		// requirements
		const float IN = state[44]; // prepared by kernel integration
		
		// exposures
		float OUT; // to be filled in
		{
			// fixed properties
			const float PARM = constants[44];
			
			// state variables
			const float STATE = state[2412];
			
			// derived variables are a function of contatnts, states, and requirements\
			const float DERIVED = PARM * STATE * IN;
			OUT = DERIVED; // is also an exposure
			
			// that's all for, if you want to update state:
			float derivative_for_STATE = DERIVED; // time derivtive for STATE
			state[2412] = state[2412] + dt * derivative_for_STATE; // or use a pointer, C has no references or lambdas :c
			
			// or initialize, or switch state:
			state[2412] = DERIVED * 3; // OnStart
			if( DERIVED < 3 ){
				state[2412] = DERIVED * 3; // OnCondition
			}
			
		}
	
	Workflow for singleton types:
		allocate
		and independently:
			- write exposure code (with assigned)
			- write update code   (with assigned)
		and
			- determine constants, emplace once on allocated signature
		
	Workflow for multitude types:
		allocate tables
		and independently:
			- write exposure code
			- write update code
		and
			- determine prototype constants, to be inherited by instances
			- extend tables by the eventual amount of instances
	
	*/
	struct SignatureAppender_Single : public ISignatureAppender{
			
		virtual size_t Constant( Real default_value, const std::string &for_what ) const {
			size_t Index = wig.constants.size();
			wig.constants.push_back(default_value);
			wig.constants_names[Index] = for_what;
			return Index;
		}
		virtual size_t StateVariable( Real default_value, const std::string &for_what ) const {
			size_t Index = wig.state.size();
			wig.state.push_back(default_value);
			wig.state_names[Index] = for_what;
			return Index;
		}
		virtual size_t Constant( std::vector<Real> bunch_of_values, const std::string &for_what ) const {
			size_t Index = wig.constants.size();
			AppendToVector(wig.constants, bunch_of_values);
			wig.constants_names[Index] = for_what;
			return Index;
		}
		virtual size_t StateVariable( std::vector<Real> bunch_of_values, const std::string &for_what ) const {
			size_t Index = wig.state.size();
			AppendToVector(wig.state, bunch_of_values);
			wig.state_names[Index] = for_what;
			return Index;
		}
		
		virtual std::string ReferTo_Const( size_t index ) const {
			return "local_constants["+itos(index)+"]";
		}
		virtual std::string ReferTo_State( size_t index ) const {
			return "local_state["+itos(index)+"]";
		}
		virtual std::string ReferTo_StateNext( size_t index ) const {
			return "local_stateNext["+itos(index)+"]";
		}
		
		CellInternalSignature::WorkItemDataSignature &wig;
		
		SignatureAppender_Single( CellInternalSignature::WorkItemDataSignature &_w )
		: wig(_w){ }
	};
	struct SignatureAppender_Table : public ISignatureAppender{
		
		virtual size_t Constant( Real default_value, const std::string &for_what ) const {
			size_t Index = wig.tables_const_f32.size();
			wig.tables_const_f32.push_back( CellInternalSignature::TableInfo(for_what) );
			
			// do something with the table's default variable
			wig.prototype_const.push_back(default_value);
			return Index;
		}
		virtual size_t StateVariable( Real default_value, const std::string &for_what ) const {
			size_t Index = wig.tables_state_f32.size();
			wig.tables_state_f32.push_back( CellInternalSignature::TableInfo(for_what) );
			// tables_for_synapse_post.push_back(table_G); they will be taken from LEMS const/state list, anyway
			
			// do something with the table's default variable
			wig.prototype_state.push_back(default_value);
			return Index;
		}
		size_t ConstI64( const std::string &for_what ) const {
			size_t Index = wig.tables_const_i64.size();
			wig.tables_const_i64.push_back( CellInternalSignature::TableInfo(for_what) );
			return Index;
		}
		size_t StateI64( const std::string &for_what ) const {
			size_t Index = wig.tables_state_i64.size();
			wig.tables_state_i64.push_back( CellInternalSignature::TableInfo(for_what) );
			return Index;
		}
		
		// but it could be there is no default valus
		size_t Constant( const std::string &for_what ) const {
			return Constant( NAN, for_what );
		}
		size_t StateVariable( const std::string &for_what ) const {
			return StateVariable( NAN, for_what );
		}
		
		virtual std::string ReferTo_Const( size_t index ) const {
			return "local_const_table_f32_arrays["+itos(index)+"][instance]";
		}
		virtual std::string ReferTo_State( size_t index ) const {
			return "local_state_table_f32_arrays["+itos(index)+"][instance]";
		}
		virtual std::string ReferTo_StateNext( size_t index ) const {
			return "local_stateNext_table_f32_arrays["+itos(index)+"][instance]";
		}
		CellInternalSignature::WorkItemDataSignature &wig;
		SignatureAppender_Table( CellInternalSignature::WorkItemDataSignature &_w )
		: wig(_w) { }
	};
	
	// Note: symbol values are assumed to be in engine native scales
	// TODO roll the random call counter into the random() functor, unless we want to keep a certain traversal order
	//  ... which is still better enforced by calling the functor
	auto EvaluateLemsExpression = []( 
		const TermTable &expression, const auto &symbol_values, const auto &symbol_dimensions, 
		const DimensionSet &dimensions, Int &random_call_counter, 
		Real &val_out, Dimension &dim_out, const auto &call_random 
	){
		
		auto Evaluate = [ ]( 
			const auto &Evaluate, // the most elegant way to make a lambda recursive, oh well
			const TermTable &expression, int node, const auto &symbol_values, const auto &symbol_dimensions, 
			const DimensionSet &dimensions, Int &random_call_counter, 
			Real &val_out, Dimension &dim_out, const auto &call_random 
		)
		-> void // Return type must be made explicit before recursive calls in the code, a line with if(0) return; could also work
		{
			
			const auto &tab = expression;
			//printf("node %d \n", node);
			auto &term = tab[node];
			if(term.type == Term::VALUE){
				val_out = term.value;
				dim_out = Dimension::Unity();
			}
			else if(term.type == Term::SYMBOL){
				val_out = symbol_values[term.symbol];
				dim_out = symbol_dimensions[term.symbol];
			}
			else if(term.isUnary() || term.isBinaryOperator()){
				
				// it is a math operation, could cause  dimension change -> possible conversion factor to shift between engine units
				LemsUnit conversion_factor =  dimensions.GetNative(Dimension::Unity()); // or override
				
				if(term.isBinaryOperator()){
					Real      val_l, val_r;
					Dimension dim_l, dim_r;
					
					Evaluate(Evaluate, expression, term.left , symbol_values, symbol_dimensions, dimensions, random_call_counter, val_l, dim_l, call_random );
					Evaluate(Evaluate, expression, term.right, symbol_values, symbol_dimensions, dimensions, random_call_counter, val_r, dim_r, call_random );
					
					if(term.type == Term::PLUS){
						val_out = val_l + val_r;
						dim_out = dim_r; // NeuroML API should have ensured consistency
					}
					else if(term.type == Term::MINUS ){
						val_out = val_l - val_r;
						dim_out = dim_r; // NeuroML API should have ensured consistency
					}
					else if(term.type == Term::TIMES ){
						val_out = val_l * val_r;
						dim_out = dim_l * dim_r;
						// correct for change of engine units...
						val_out = ( dimensions.GetNative(dim_l) * dimensions.GetNative(dim_r) ).ConvertTo(val_out, dimensions.GetNative(dim_out)) ;
					}
					else if(term.type == Term::DIVIDE){
						val_out = val_l / val_r;
						dim_out = dim_l / dim_r;
						// correct for change of engine units...
						val_out = ( dimensions.GetNative(dim_l) / dimensions.GetNative(dim_r) ).ConvertTo(val_out, dimensions.GetNative(dim_out)) ;
					}
					else if(term.type == Term::LT    ){
						val_out = ( val_l < val_r ) ? 1 : 0;
						dim_out = Dimension::Unity(); // boolean
					}
					else if(term.type == Term::GT    ){
						val_out = ( val_l > val_r ) ? 1 : 0;
						dim_out = Dimension::Unity(); // boolean
					}
					else if(term.type == Term::LEQ   ){
						val_out = ( val_l <= val_r ) ? 1 : 0;
						dim_out = Dimension::Unity(); // boolean
					}
					else if(term.type == Term::GEQ   ){
						val_out = ( val_l >= val_r ) ? 1 : 0;
						dim_out = Dimension::Unity(); // boolean
					}
					else if(term.type == Term::EQ    ){
						val_out = ( val_l == val_r ) ? 1 : 0;
						dim_out = Dimension::Unity(); // boolean
					}
					else if(term.type == Term::NEQ   ){
						val_out = ( val_l != val_r ) ? 1 : 0;
						dim_out = Dimension::Unity(); // boolean
					}
					else if(term.type == Term::AND   ){
						val_out = ( (bool)val_l && (bool)val_r ) ? 1 : 0;
						dim_out = Dimension::Unity(); // boolean
					}
					else if(term.type == Term::OR    ){
						val_out = ( (bool)val_l || (bool)val_r ) ? 1 : 0;
						dim_out = Dimension::Unity(); // boolean
					}
					else if(term.type == Term::POWER ){
						val_out = std::pow( val_l, val_r );
						dim_out = Dimension::Unity(); // XXX i know
					}
					else{
						assert(false);
					}
				}
				else if( term.type == Term::RANDOM ){
					Real val_r;
					Dimension dim_r;
					Evaluate(Evaluate, expression, term.right, symbol_values, symbol_dimensions, dimensions, random_call_counter, val_r, dim_r, call_random );
					
					val_out = call_random(val_r, random_call_counter);
					random_call_counter++;
					
					dim_out = dim_r; // should be pure number
				}
				else if(term.isUnaryFunction()){
					Real val_r;
					Dimension dim_r;
					Evaluate(Evaluate, expression, term.right, symbol_values, symbol_dimensions, dimensions, random_call_counter, val_r, dim_r, call_random );
					
							if(term.type == Term::ABS   ){ val_out = std::abs  ( val_r ); }
					else if(term.type == Term::SQRT  ){ val_out = std::sqrt ( val_r ); }
					else if(term.type == Term::SIN   ){ val_out = std::sin  ( val_r ); }
					else if(term.type == Term::COS   ){ val_out = std::cos  ( val_r ); }
					else if(term.type == Term::TAN   ){ val_out = std::tan  ( val_r ); }
					else if(term.type == Term::SINH  ){ val_out = std::sinh ( val_r ); }
					else if(term.type == Term::COSH  ){ val_out = std::cosh ( val_r ); }
					else if(term.type == Term::TANH  ){ val_out = std::tanh ( val_r ); }
					else if(term.type == Term::EXP   ){ val_out = std::exp  ( val_r ); }
					else if(term.type == Term::LOG10 ){ val_out = std::log10( val_r ); }
					else if(term.type == Term::LN    ){ val_out = std::log  ( val_r ); }
					else if(term.type == Term::CEIL  ){ val_out = std::ceil ( val_r ); }
					else if(term.type == Term::FLOOR ){ val_out = std::floor( val_r ); }
					else if(term.type == Term::HFUNC ){ val_out = ( val_r < 0 ) ? 0 : 1; }
					else{
						assert(false);
					}
					
					// math functions conveniently map from pure number to pure number, LATER specialcase if needed
					dim_out = dim_r;
				}
				else if(term.isUnaryOperator()){
					Real val_r;
					Dimension dim_r;
					Evaluate(Evaluate, expression, term.right, symbol_values, symbol_dimensions, dimensions, random_call_counter, val_r, dim_r, call_random );
					
					if(term.type == Term::UMINUS ){
						val_out = -val_r;
						dim_out = dim_r;
					}
					else if(term.type == Term::UPLUS ){
						val_out = +val_r;
						dim_out = dim_r;
					}
					else if(term.type == Term::NOT ){
						val_out = (bool)val_r;
						dim_out = Dimension::Unity(); // boolean
					}
					else{
						assert(false);
					}
				}
				else{
					assert(false);
				}
				
			}
			else{
				printf("unknown term %d !\n", term.type);
				assert(false); // TODO something better
			}
			// append debug info for e.g. dimensions
			// printf("%g %s\n", val_out, dimensions.Stringify(dim_out).c_str());
			//printf("bye node %d \n", node);
		};
		//TermTable::printTree(expression.tab, expression.tab.expression_root, 0);
		val_out = NAN; // just to initialize
		dim_out = Dimension::Unity(); // just to initialize
		Evaluate(Evaluate, expression, expression.expression_root, symbol_values, symbol_dimensions, dimensions, random_call_counter, val_out, dim_out, call_random);
		return;
	};
	
	struct DescribeLems{
		
		static std::string ExpressionInfix( const ComponentType::ResolvedTermTable &expression, const ComponentType &type, const DimensionSet &dimensions, Int &random_call_counter, Dimension &dim_out ){
			struct Help{
				static void Infix(
					const ComponentType::ResolvedTermTable &expression, int node,
					const ComponentType type, const DimensionSet &dimensions, Int &random_call_counter,
					std::string &out, Dimension &dim_out
				){
					const auto &tab = expression.tab;
					//printf("node %d \n", node);
					auto &term = tab[node];
					if(term.type == Term::VALUE){
						out += "( (float)"+accurate_string(term.value)+")";
						dim_out = Dimension::Unity();
					}
					else if(term.type == Term::SYMBOL){
						auto assigned_seq = expression.resolved[term.symbol];
						out += "*Lems_assigned_";
						out += itos(assigned_seq);
						dim_out = type.getNamespaceEntryDimension(assigned_seq);
					}
					else if(term.isUnary() || term.isBinaryOperator()){
						
						// it is a math operation, could cause  dimension change -> possible conversion factor to shift between engine units
						LemsUnit conversion_factor =  dimensions.GetNative(Dimension::Unity()); // or override
						out += "( "; // for conversion factor
						const static std::map< Term::Type, const char * > term_strings = {
							
							{ Term::LEQ   , "<=" },
							{ Term::GEQ   , ">=" },
							{ Term::LT    , "<"  },
							{ Term::GT    , ">"  },
							{ Term::EQ    , "==" },
							{ Term::NEQ   , "!=" },
							{ Term::AND   , "&&" },
							{ Term::OR    , "||" },
							
							{ Term::UMINUS, "-" },
							{ Term::UPLUS , "+" },
							{ Term::NOT   , "!" },
							
							{ Term::PLUS  , "+" },
							{ Term::MINUS , "-" },
							{ Term::TIMES , "*" },
							{ Term::DIVIDE, "/" },
							{ Term::POWER ,  "powf" },
							
							{ Term::ABS   ,     "fabs" },
							{ Term::SQRT  , "sqrtf" },
							{ Term::SIN   ,  "sinf" },
							{ Term::COS   ,  "cosf" },
							{ Term::TAN   ,  "tanf" },
							{ Term::SINH  , "sinhf" },
							{ Term::COSH  , "coshf" },
							{ Term::TANH  , "tanhf" },
							{ Term::EXP   ,  "expf" },
							{ Term::LOG10 ,"log10f" },
							{ Term::LN    ,  "logf" },
							{ Term::CEIL  , "ceilf" },
							{ Term::FLOOR ,"floorf" },
							{ Term::RANDOM,"randof" },
							{ Term::HFUNC , "stepf" },
						};
						
						auto it = term_strings.find(term.type);
						if( it == term_strings.end() ){
							printf("unknown termstring  !\n");
							//assert(false);
						}
						const auto &termstr = it->second;
						if(term.isBinaryOperator()){
							Dimension dim_l, dim_r;
							if(term.type == Term::POWER){
								out += termstr;
								out += "( ";
								Infix(expression, term.left, type, dimensions, random_call_counter, out, dim_l);
								out += " , ";
								Infix(expression, term.right, type, dimensions, random_call_counter, out, dim_r);
								out += " )";
							}
							else{
								// generally, use infix
								out += "( ";
								Infix(expression, term.left, type, dimensions, random_call_counter, out, dim_l);
								out += " ";
								out += termstr;
								out += " ";
								Infix(expression, term.right, type, dimensions, random_call_counter, out, dim_r);
								out += " )";
							}
							if(term.type == Term::POWER){
								dim_out = Dimension::Unity(); // TODO plead to the community to come along in banning usage of pow sqrt and such with non pure numbers - otherwise the dimensionality calculus becomes extremely complicated and yet easy to break. How does Brian handle it anyway?
							}
							else if(term.type == Term::PLUS){
								dim_out = dim_r; // NeuroML API should have ensured consistency
							}
							else if(term.type == Term::MINUS){
								dim_out = dim_r; // NeuroML API should have ensured consistency
							}
							else if(term.type == Term::TIMES){
								dim_out = dim_l * dim_r;
								conversion_factor = ( dimensions.GetNative(dim_l) * dimensions.GetNative(dim_r) ).to( dimensions.GetNative(dim_out) );
							}
							else if(term.type == Term::DIVIDE){
								dim_out = dim_l / dim_r;
								conversion_factor = ( dimensions.GetNative(dim_l) / dimensions.GetNative(dim_r) ).to( dimensions.GetNative(dim_out)) ;
							}
							
						}
						else if( term.type == Term::RANDOM ){
							out += termstr;
							out +="( ";
							Dimension dim_r;
							Infix(expression, term.right, type, dimensions, random_call_counter, out, dim_r);
							out += " , rng_object_id, instance, step, rng_offset + "; // TODO use a macro, cleaner this way
							out += itos( random_call_counter );
							random_call_counter++;
							out += " )";
							dim_out = dim_r; // should be pure number
						}
						else if(term.isUnaryFunction()){
							
							out += termstr;
							out += "( ";
							Dimension dim_r;
							Infix(expression, term.right, type, dimensions, random_call_counter, out, dim_r);
							out += " )";
							
							// math functions conveniently map from pure number to pure number, LATER specialcase if needed
							dim_out = dim_r;
						}
						else if(term.isUnaryOperator()){
							auto termstr = it->second;
							out += "( ";
							out += termstr;
							out += " ";
							Dimension dim_r;
							Infix(expression, term.right, type, dimensions, random_call_counter, out, dim_r);
							out += " )";
							
							// unary operators map from and to same dimension, LATER specialcase if needed
							dim_out = dim_r;
						}
						else{
							// who knows
							out += " ??? ";
							dim_out = Dimension::Unity();
						}
						
						// apply scale adjustments for engine units
						// only multiplication, division lead to unit conversion, for now
						// could possibly have something like a unary square operator etc.
						// though it really is best to keep a consistent units system, so factors are not required
						// otherwise what should be done if no native unit is specified for an intermediate result ???
						// what unit should be two differently-occurring addends be normalized at?
						// -> get a fallback basis of fundamental units, to avoid extreme rescaling TODO
						out += Convert::Suffix(conversion_factor);
						out += " )"; // for conversion factor
					}
					else{
						printf("unknown term !\n");
						assert(false); // TODO something better
					}
					// append debug info for e.g. dimensions
					out += "/* ";
					out += dimensions.Stringify(dim_out);
					out += " */";
					//printf("bye node %d \n", node);
				}
			};
			std::string ret;
			//TermTable::printTree(expression.tab, expression.tab.expression_root, 0);
			dim_out = Dimension::Unity(); // just to initialize
			Help::Infix(expression, expression.tab.expression_root, type, dimensions, random_call_counter, ret, dim_out);
			return ret;
		}
		static std::string ExpressionInfix( const ComponentType::ResolvedTermTable &expression, const ComponentType &type, const DimensionSet &dimensions, Int &random_call_counter){
			Dimension dim_out;
			return ExpressionInfix(expression, type, dimensions, random_call_counter, dim_out);
		}
		static std::string GetExposureVar(const ComponentType &type, Int exp_seq){
			std::string ret = "Lems_";
			const auto &exp = type.exposures.get(exp_seq);
			if(exp.type == ComponentType::Exposure::STATE){
				ret += "state";
			}
			else if(exp.type == ComponentType::Exposure::DERIVED){
				ret += "derived";
			}
			else{
				ret += "unknown";
			}
			ret += "_";
			ret += itos(exp.seq);
			
			return ret;
		};
		static CellInternalSignature::ComponentValueInstance GetValues(const ComponentType &type, const ComponentInstance &instance){
			CellInternalSignature::ComponentValueInstance ret;
			
			std::vector<Real> customized_constants(type.properties.contents.size());
			// try overriding the properties with some parms for this component instance
			for(size_t seq = 0; seq < type.properties.contents.size(); seq++ ){
				customized_constants[seq] = type.properties.get( (Int)seq ).value;
			}
			for(auto parm : instance.parms){
				customized_constants[parm.seq] = parm.value;
			}
			
			// fill initial state with NANs, let init steps sort the OnStart assignments and derived state variables out
			// actually fill them with zero, because that's the often relied-upon undocumented behaviour
			std::vector<Real> customized_initstates(type.state_variables.contents.size(), 0);
			
			ret.properties = customized_constants;
			ret.statevars = customized_initstates;
			
			return ret;
		}
		static CellInternalSignature::ComponentSubSignature AllocateSignature(const ComponentType &type, const ComponentInstance &instance, const ISignatureAppender *Add, std::string for_what){
			// TODO random variables should be examined here, for consideration MUCH LATER
			CellInternalSignature::ComponentSubSignature ret;
			// TODO rework to split values from allocation
			
			auto vals = GetValues( type, instance );
			
			// the persistent numbers are parameters(of any type), and state variables
			for(size_t seq = 0; seq < type.properties.contents.size(); seq++ ){
				size_t Index = Add->Constant( vals.properties[seq], for_what + std::string(" Property ") + itos(seq) ); // TODO dimensions + " ("++")"
				
				ret.properties_to_constants.push_back(CellInternalSignature::ComponentSubSignature::Entry(Index, CellInternalSignature::ComponentSubSignature::Entry::ValueType::F32));
			}
			for(size_t seq = 0; seq < type.state_variables.contents.size(); seq++ ){
				size_t Index = Add->StateVariable( vals.statevars[seq], for_what + std::string(" State ") + itos(seq) );
				
				ret.statevars_to_states.push_back({Index, CellInternalSignature::ComponentSubSignature::Entry::ValueType::F32});
			}
			
			return ret;
		}
		/*
		
		static void ApplySubsigToInstance(const ComponentType &type, const CellInternalSignature::ComponentSubSignature &subsig, const ISignatureAppender *Add){
			// mutatesparametrization of components; singletons and multitudes are separately synapses and inputs, other parms are semi-standard
		}
		*/
		// Assume inputs have already been defined, this is what the component assigns itself (derived etc.)
		static std::string Assigned(const ComponentType &type, const DimensionSet &dimensions, const CellInternalSignature::ComponentSubSignature &subsig, const ISignatureAppender *Add, const std::string &for_what, const std::string &line_prefix, Int &random_call_counter, bool debug = false){
			const auto &tab = line_prefix; // for a more convenient name
			char tmps[2000];
			std::string ret;
			
			// TODO replace common exposures by LEMS names using a persistent name table
			const static NameMap< Int ComponentType::CommonRequirements::* > common_requirement_names = {
				{"time"					, &ComponentType::CommonRequirements::time },
				{"temperature"			, &ComponentType::CommonRequirements::temperature },
				{"Vcomp"				, &ComponentType::CommonRequirements::membrane_voltage },
				{"Acomp"				, &ComponentType::CommonRequirements::membrane_surface_area },
				{"Iion"         		, &ComponentType::CommonRequirements::ion_current },					
				{"InitConcIntra"		, &ComponentType::CommonRequirements::initial_concentration_intra },
				{"InitConcExtra"	    , &ComponentType::CommonRequirements::initial_concentration_extra },
				{"Ca_concentration"		, &ComponentType::CommonRequirements::calcium_concentration_intra },
				{"alpha"				, &ComponentType::CommonRequirements::gate_rate_alpha },
				{"beta "				, &ComponentType::CommonRequirements::gate_rate_beta  },
				{"rateScale "			, &ComponentType::CommonRequirements::gate_rate_scale },
				{"Vpeer"				, &ComponentType::CommonRequirements::peer_voltage },
				{"block_factor"			, &ComponentType::CommonRequirements::block_factor },
				{"plasticity_factor"	, &ComponentType::CommonRequirements::plasticity_factor },
				{"plasticity_factor"	, &ComponentType::CommonRequirements::plasticity_factor },
				{"external_current"		, &ComponentType::CommonRequirements::external_current },
			};
			// sort the lines, to protect the maintainer's psyche
			std::vector<std::string> req_lines;
			for( auto keyval : common_requirement_names ){
				auto name = keyval.first;
				auto ptr = keyval.second;
				Int req_seq = type.common_requirements.*ptr;
				if(req_seq >= 0){
					sprintf(tmps, "float Lems_requirement_%ld = %s;\n", req_seq, name); 
					//ret += tab+tmps;
					req_lines.push_back( tab+tmps );
				}
			}
			std::sort(req_lines.begin(), req_lines.end());
			for(auto line : req_lines) ret += line;
			
			// perhaps split off to different loop LATER
			const static NameMap< Int ComponentType::CommonEventInputs::* > common_eventin_names = {
				{"spike_in_flag"			, &ComponentType::CommonEventInputs::spike_in },
			};
			std::vector<std::string> eventin_lines;
			for( auto keyval : common_eventin_names ){
				auto name = keyval.first;
				auto ptr = keyval.second;
				Int req_seq = type.common_event_inputs.*ptr;
				if(req_seq >= 0){
					sprintf(tmps, "char Lems_eventin_%ld = %s;\n", req_seq, name); 
					//ret += tab+tmps;
					eventin_lines.push_back( tab+tmps );
				}
			}
			std::sort(eventin_lines.begin(), eventin_lines.end());
			for(auto line : eventin_lines) ret += line;
			
			// set vars for event output flags
			for(size_t i = 0; i < type.event_outputs.contents.size(); i++){
				sprintf(tmps, "float Lems_evout_%zd = 0;", i );
				ret += tab+tmps+"\n";
			}
			
			
			// TODO clarify better by keeping the actual names
			ret += tab+"// fixed properties "+for_what+"\n";
			for(size_t i = 0; i < type.properties.contents.size(); i++){
				sprintf(tmps, "float Lems_property_%zd = %s;", i, Add->ReferTo_Const(subsig.properties_to_constants.at(i).index).c_str() );
				ret += tab+tmps+"\n";
			}
			ret += tab+"// state variables "+for_what+"\n";
			for(size_t i = 0; i < type.state_variables.contents.size(); i++){
				sprintf(tmps, "float Lems_state_%zd = %s;", i, Add->ReferTo_State(subsig.statevars_to_states.at(i).index).c_str() );
				ret += tab+tmps+"\n";
			}
			
			ret += tab+"// declare derived variables "+for_what+"\n";
			for(size_t i = 0; i < type.derived_variables.contents.size(); i++){
				sprintf(tmps, "float Lems_derived_%zd = NAN;", i);
				ret += tab+tmps+"\n";
			}
			
			ret += tab+"// common read-only namespace? "+for_what+"\n";
			// TODO eliminate, for less confusion to readers and compilers
			for(size_t i = 0; i < type.name_space.contents.size(); i++){
				sprintf(tmps, "float *Lems_assigned_%zd = &Lems_", i); 
				ret += tab+tmps;
				auto namet = type.name_space.get(i).type;
				if(namet == ComponentType::NamespaceThing::PROPERTY){
					ret += "property";
				}
				else if(namet == ComponentType::NamespaceThing::REQUIREMENT){
					ret += "requirement";
				}
				else if(namet == ComponentType::NamespaceThing::STATE){
					ret += "state";
				}
				else if(namet == ComponentType::NamespaceThing::DERIVED){
					ret += "derived";
				}
				
				sprintf(tmps, "_%ld;\n", type.name_space.get(i).seq);
				ret += tmps;
			}
			
			ret += tab+"// compute derived "+for_what+"\n";
			for(size_t i = 0; i < type.derived_variables_topological_order.size(); i++){
				Int seq = type.derived_variables_topological_order[i];
				const auto &dervar = type.derived_variables.get(seq);
				
				if(dervar.type == ComponentType::DerivedVariable::VALUE){
					sprintf(tmps, "Lems_derived_%ld = ", seq);
					assert( dervar.cases.size() == 0 );
					auto expression_string = ExpressionInfix(dervar.value, type, dimensions, random_call_counter);
					ret += tab+tmps+expression_string+";\n";
				}
				else if(dervar.type == ComponentType::DerivedVariable::CONDITIONAL){
					sprintf(tmps, "Lems_derived_%ld = 0;", seq); ret += tab + tmps;
					
					ret += tab + "if( 0 );\n"; // to avoid extra logic for the first 'if' and the case only 'default' case exists
					// conditional cases
					for( size_t case_seq = 0; case_seq < dervar.cases.size(); case_seq++ ){
						const auto &deri_case = dervar.cases[case_seq];
						if( (Int)case_seq == dervar.default_case ) continue;
						
						auto condition_string = ExpressionInfix(deri_case.condition, type, dimensions, random_call_counter);
						ret += tab + "else if( " + condition_string + " ){\n";
						
						auto value_string = ExpressionInfix(deri_case.value, type, dimensions, random_call_counter);
						sprintf(tmps, "\tLems_derived_%ld = ", seq);
						ret += tab + tmps + value_string + ";\n";
						
						ret += tab + "}\n";
					}
					// and default case
					if( dervar.default_case >= 0 ){
						const auto &deri_case = dervar.cases[dervar.default_case];
						
						ret += tab + "else{\n";
						
						auto value_string = ExpressionInfix(deri_case.value, type, dimensions, random_call_counter);
						sprintf(tmps, "\tLems_derived_%ld = ", seq);
						ret += tab + tmps + value_string + ";\n";
						
						ret += tab + "}\n";
					}
				}
				else{
					sprintf(tmps, "internal error: assigned derived variable %ld type %d\n", seq, dervar.type);
					return tmps; // why not
				}
				// if(debug){
				// 	sprintf(tmps, "printf(\"derived %zd = %%.17g \\n\", Lems_derived_%zd);\n", seq, seq); ret += tab+tmps;
				// }
			}
			
			return ret;
		}
		// NOTE exposures should be added after dynamics, because spikes outputs are triggered on OnCondition atatements !
		// If there's a need for getting them before dynamics, add a second scan for OnConditions into Assigned section LATER
		static std::string Exposures(const ComponentType &type, const std::string &for_what, const std::string &line_prefix, bool debug = false){
			const auto &tab = line_prefix; // for a more convenient name
			//char tmps[2000];
			std::string ret;
			ret += tab+"// exposures "+for_what+"\n";
			
			std::vector<std::string> exposure_lines;
			// XXX sort the keys first, because hashtable is unordered; LATER if it occurs or matters I guess?
			for( auto keyval : ComponentType::CommonExposures::names ){
				auto name = keyval.first;
				auto ptr = keyval.second;
				if(type.common_exposures.*ptr >= 0){
					exposure_lines.push_back( "float Lems_exposure_"+std::string(name)+" = " + GetExposureVar(type, type.common_exposures.*ptr)+";\n" );
				}
			}
			std::sort(exposure_lines.begin(), exposure_lines.end());
			for(auto line : exposure_lines) ret += tab+line;
			
			std::vector<std::string> eventout_lines;
			for( auto keyval : ComponentType::CommonEventOutputs::names ){
				auto name = keyval.first;
				auto ptr = keyval.second;
				if(type.common_event_outputs.*ptr >= 0){
					eventout_lines.push_back( "char Lems_eventout_"+std::string(name)+" = Lems_evout_"+itos(type.common_event_outputs.*ptr)+";\n" );
				}
			}
			std::sort(eventout_lines.begin(), eventout_lines.end());
			for(auto line : eventout_lines) ret += tab+line;
			
			return ret;
		}
		// Assume assigned values have already been defined, this updates state variables (rates, conditions etc.)
		static std::string Update(const ComponentType &type, const DimensionSet &dimensions, const CellInternalSignature::ComponentSubSignature &subsig, const ISignatureAppender *Add, const std::string &for_what, const std::string &line_prefix, Int &random_call_counter, bool debug = false){
			
			const auto &tab = line_prefix; // for a more convenient name
			char tmps[2000];
			std::string ret;
			
			// TODO A sequence of AssignState commands may have state variables assigned with values of other, previously assigned state variables (possibly even modifying derived variables !!!)
			// So subsequent AssignState commands should keep track of which state variables were modified, and use the updated values.
			// On the other hand, no other NeuroML implementation seems to correct derived variables for the updated value, so there's that.
			std::vector< std::vector<int> > statevar_to_assigned( type.state_variables.contents.size() );
			for(int i = 0; i < (int)type.name_space.contents.size(); i++){
				auto namet = type.name_space.get(i).type;
				auto ref_seq = type.name_space.get(i).seq;
				if(namet == ComponentType::NamespaceThing::PROPERTY){
					// skip
				}
				else if(namet == ComponentType::NamespaceThing::REQUIREMENT){
					// skip
				}
				else if(namet == ComponentType::NamespaceThing::STATE){
					statevar_to_assigned[ref_seq].push_back(i);
				}
				else if(namet == ComponentType::NamespaceThing::DERIVED){
					// TODO
				}
				
			}
			
			auto Emit_AssignState = [&](const ComponentType::StateAssignment &assign){
				auto state_seq = assign.state_seq;
				auto Index = subsig.statevars_to_states.at(state_seq).index;
				
				sprintf(tmps, "		%s = ", Add->ReferTo_StateNext(Index).c_str() );
				auto expression_string = ExpressionInfix(assign.value, type, dimensions, random_call_counter);
				ret += tab+tmps+expression_string+";\n";
				for( auto assigned_seq : statevar_to_assigned[state_seq] ){
					sprintf(tmps, "		Lems_assigned_%d = &(%s) ", assigned_seq, Add->ReferTo_StateNext(Index).c_str() );
					ret += tab+tmps+ ";\n";
				}
			};
			
			auto Emit_EventOut = [&](const ComponentType::EventOut &evout){				
				sprintf(tmps, "		Lems_evout_%ld = 1", evout.port_seq );
				ret += tab+tmps+";\n";
			};
			
			ret += tab+"if(initial_state){"+"\n";
			ret += tab+"	// initialization"+"\n";
			
			for(auto assign : type.on_start){
				Emit_AssignState(assign);
			}
			
			ret += tab+"}else{"+"\n";
			
			ret += tab+"	// dynamics"+"\n";
			ret += tab+"	// (highest up is lowest priority)"+"\n";
			
			ret += tab+"	// time derivatives"+"\n";
			for(size_t seq = 0; seq < type.state_variables.contents.size(); seq++){
				const auto &state_variable = type.state_variables.get(seq);
				auto Index = subsig.statevars_to_states.at(seq).index;
				
				if( state_variable.dynamics == ComponentType::StateVariable::DYNAMICS_NONE ){
					// keep it unchanged
					sprintf(tmps, "	%s = %s;", Add->ReferTo_StateNext(Index).c_str(), Add->ReferTo_State(Index).c_str() );
					ret += tab+tmps+"\n";
				}
				else if( state_variable.dynamics == ComponentType::StateVariable::DYNAMICS_CONTINUOUS ){
					// Euler integration, add fancier techniques LATER
					Dimension dim_of_derivative;
					
					auto expression_string = ExpressionInfix(state_variable.derivative, type, dimensions, random_call_counter, dim_of_derivative);
					
					LemsUnit conversion_factor = ( dimensions.GetNative(dim_of_derivative) * dimensions.GetNative(LEMS_Time) ).to( dimensions.GetNative(state_variable.dimension)) ;
					
					sprintf(tmps, "	float Lems_derivative_%zd = ", seq );
					ret += tab + tmps + expression_string + ";\n";
					
					
					sprintf(tmps, "	%s = %s + dt * Lems_derivative_%zd%s;", Add->ReferTo_StateNext(Index).c_str(), Add->ReferTo_State(Index).c_str(), seq, Convert::Suffix(conversion_factor).c_str() );
					ret += tab+tmps+"\n";
				}
				else{
					ret += tab+"	missing dynamics for variable "+itos(seq)+"\n";
				}
			}
			// todo run conditions AND integrate in the same step, somehow;
			// though this will require two invocations of assigned variables:
			//   one to evaluate triggers to handle,
			//   and another one to re-evaluate assigned variables AFTER triggers have run
			// NEURON does this by handling WATCH statements as NetCon events, asynchronously with integration
			// but can this be done in a less hairy manner? perhaps LATER
			// Keep in mind that NEURON's WATCH statements are rather cumbersome to use
			
			auto HandleDoStuff = [&]( const auto &oncase){
				for( auto assign : oncase.assign ){
					Emit_AssignState(assign);
				}
				
				for( auto evout : oncase.event_out ){
					Emit_EventOut(evout);
				}
			};
			
			// ret += tab+"// conditional updates, both for initialization and simulation"+"\n"; not now really, if case arises LATER
			ret += tab+"// conditional updates, during simulation"+"\n";
			for(size_t seq = 0; seq < type.on_conditions.size(); seq++){
				const auto &onco = type.on_conditions.at(seq);
				auto expression_string = ExpressionInfix(onco.test, type, dimensions, random_call_counter);
				ret += tab+"if( "+expression_string+" ){\n";
				HandleDoStuff(onco);
				ret += tab+"}"+"\n";
			}
			// perhaps these should be reordered?
			
			// perhaps split off to different loop LATER
			for(size_t seq = 0; seq < type.on_events.size(); seq++){
				const auto &onen = type.on_events.at(seq);
				sprintf(tmps, "Lems_eventin_%ld", onen.in_port_seq);
				std::string expression_string = tmps;
				ret += tab+"if( "+expression_string+" ){\n";
				HandleDoStuff(onen);
				ret += tab+"}"+"\n";
			}
			
			// TODO on inbound events (perhaps these should be followed along with updates, but what about dependencies then?)
			
			ret += tab+"}"+"\n";
			
			return ret;
		}
	};
	
	auto DescribeLems_AppendTableEntry = [ &tabs, &component_types ]( size_t work_unit, const ComponentInstance &comp_instance, const CellInternalSignature::ComponentSubSignature &subsig ){
		
		// TODO be a little dirty until it stabilizes
		// grab any constants from the property type, then replace them with instance parms, then append them to the vectors
		
		const auto off_cf32 = tabs.global_table_const_f32_index[work_unit];
		const auto off_sf32 = tabs.global_table_state_f32_index[work_unit];
		auto &tab_cf32 = tabs.global_tables_const_f32_arrays;
		auto &tab_sf32 = tabs.global_tables_state_f32_arrays;
		
		const ComponentType &comp_type = component_types.get(comp_instance.id_seq);
		
		const auto vals = DescribeLems::GetValues( comp_type, comp_instance );
		
		for( size_t seq = 0; seq < vals.properties.size(); seq++ ){
			tab_cf32[off_cf32 + subsig.properties_to_constants[seq].index].push_back(vals.properties.at(seq));
		}
		
		for( size_t seq = 0; seq < vals.statevars.size(); seq++ ){
			tab_sf32[off_sf32 + subsig.statevars_to_states[seq].index].push_back(vals.statevars.at(seq));
		}
		
	};
	
	
	struct InlineLems_AllocatorCoder{
		const Model &model;
		Int &random_call_counter; // XXX convert to allocator
		const SignatureAppender_Single &AppendSingle;
		const SignatureAppender_Table &AppendMulti;
		
		// TODO return bool for error handling
		
		std::string SingleInstance(
			const ComponentInstance &compinst,
			const std::string &tab, const std::string &for_what,
			CellInternalSignature::ComponentSubSignature &component, bool debug = false
		) const {
			std::string code;
			const auto &comptype = model.component_types.get(compinst.id_seq);
			
			component = DescribeLems::AllocateSignature(comptype, compinst, &AppendSingle, for_what + " LEMS");
			
			
			code += tab+"// LEMS component\n";
			std::string lemscode = DescribeLems::Assigned(comptype, model.dimensions, component, &AppendSingle, for_what, tab, random_call_counter, debug );
			code += lemscode;
			
			// also add integration code here, to finish with component code (and get event outputs !)
			code += tab+"// integrate inline\n";
			std::string lemsupdate = DescribeLems::Update(comptype, model.dimensions, component, &AppendSingle, for_what, tab, random_call_counter, debug);
			code += lemsupdate;
			
			code += tab+"// expose inline\n";
			code += DescribeLems::Exposures( comptype, for_what, tab, debug );
			
			return code;
		}
		
		std::string TableInstances(
			const std::string &tab, const std::string &for_what,
			CellInternalSignature::ComponentSubSignature &compsubsig
		) const {
			std::string code;
		
			// if no constant(ie variable property) tables exist for this LEMS type(could be!), 
			//  add a dummy table if no constant tables exist, just to keep track of no.of instances
			// optimize that LATER when instance pop.counters are added as separate entities
			if( compsubsig.properties_to_constants.empty() ){
				// make this prettier LATER
				std::size_t Index = (&AppendMulti)->Constant( NAN, for_what + std::string(" Dummy Property") );
				
				compsubsig.properties_to_constants.push_back(CellInternalSignature::ComponentSubSignature::Entry(Index, CellInternalSignature::ComponentSubSignature::Entry::ValueType::F32));
			}
			code += tab+"const long long Instances = local_const_table_f32_sizes["
					+ itos(compsubsig.properties_to_constants[0].index)
					+ "]; //same for all parallel arrays\n";
					
			return code;
		}
		
		std::string TableLoop(
			const std::string &tab, const std::string &for_what,
			CellInternalSignature::ComponentSubSignature &compsubsig // LATER const
		) const {
			std::string code;
			
			code += TableInstances( tab, for_what, compsubsig );
			code += tab+"for(long long instance = 0; instance < Instances; instance++)\n";
			
			return code;
		}
		
		std::string TableInner(
			const std::string &tab, const std::string &for_what,
			const ComponentType &comptype,
			const CellInternalSignature::ComponentSubSignature &compsubsig,
			const std::string &requirement_code, const std::string &exposure_code, bool debug = false
		) const {
			std::string code;
			{
			const auto &bat = tab;
			std::string tab = bat + "\t";
			
			code += tab+"// External Requirements\n";
			code += tab+requirement_code;
			
			code += tab+"// LEMS component\n";
			std::string lemscode = DescribeLems::Assigned(comptype, model.dimensions, compsubsig, &AppendMulti, for_what, tab, random_call_counter, debug );
			code += lemscode;
			
			code += tab+"// integrate inline\n";
			std::string lemsupdate = DescribeLems::Update(comptype, model.dimensions, compsubsig, &AppendMulti ,for_what, tab, random_call_counter, debug);
			code += lemsupdate;
			}
			
			code += tab+"// expose inline\n";
			code += DescribeLems::Exposures( comptype, for_what, tab, debug );
			code += tab+"// External Exposures\n";
			code += tab+exposure_code;
			
			return code;
		}
		
		std::string TableFull(
			const ComponentInstance &compinst,
			const std::string &tab, const std::string &for_what,
			CellInternalSignature::ComponentSubSignature &compsubsig,
			const std::string &requirement_code, const std::string &exposure_code, bool debug = false
		) const {
			std::string code;
			// TODO validate id_seq, to de-duplicate checking code
			const auto &comptype = model.component_types.get(compinst.id_seq);
			
			compsubsig = DescribeLems::AllocateSignature(comptype, compinst, &AppendMulti, for_what + " LEMS");
			
			code += TableLoop( tab, for_what, compsubsig );
			code += tab+"{\n";
				code += TableInner( tab, for_what, comptype, compsubsig, requirement_code, exposure_code, debug );
			code += tab+"}\n";
			
			return code;
		}
		
		InlineLems_AllocatorCoder( const Model &_m, Int &_cc, const SignatureAppender_Single &_as, const SignatureAppender_Table &_am )
			: model(_m), random_call_counter(_cc), AppendSingle(_as), AppendMulti(_am) {
				
		}
	};
	
	// kernel file boilerplate code
		
	auto EmitKernelFileHeader = [ &config ]( std::string &code ){
		(void) config; // just in case
		
		char tmps[1000];
		code += "// Generated code block BEGIN\n";
		// code += "#include <stdatomic.h>\n";
		code += "#define M_PI       3.14159265358979323846\n";
		code += "#include <math.h>\n";
		if(config.debug){
			code += "#include <stdio.h>\n";
		}
		sprintf(tmps, "typedef float * __restrict__ __attribute__((align_value (%zd))) Table_F32;\n", RawTables::ALIGNMENT); code += tmps;
		//TODO fix const correctness !
		sprintf(tmps, "typedef long long * __restrict__ __attribute__((align_value (%zd))) Table_I64;\n", RawTables::ALIGNMENT); code += tmps;
		
		// 32bit float <-> int type smuggling
		code += "typedef union { int i32; float f32; } TypePun_I32F32; typedef char static_assert[ sizeof(int) == sizeof(float) ];\n";
		code += "static float EncodeI32ToF32( int   i ){ TypePun_I32F32 cast; cast.i32 = i; return cast.f32;}\n";
		code += "static int   EncodeF32ToI32( float f ){ TypePun_I32F32 cast; cast.f32 = f; return cast.i32;}\n";
		
		// the humble but mighty step function
		code += "static float stepf( float x ){ if( x < 0 ) return 0; else return 1;  }\n";
		
		// Hash-based RNG
		code += "\n";		
		code += "// Credits to Thomas T. Wang: wang@cup.hp.com\n";
		code += "static unsigned long long hash64shift( unsigned long long key ){\n";
		code += "	key = (~key) + (key << 21); // key = (key << 21) - key - 1;\n";
		code += "	key = key ^ (key >> 24);\n";
		code += "	key = (key + (key << 3)) + (key << 8); // key * 265\n";
		code += "	key = key ^ (key >> 14);\n";
		code += "	key = (key + (key << 2)) + (key << 4); // key * 21\n";
		code += "	key = key ^ (key >> 28);\n";
		code += "	key = key + (key << 31);\n";
		code += "	return key;\n";
		code += "}\n";
		code += "static unsigned long long hash_128_to_64( unsigned long long hi, unsigned long long lo ){\n";
		code += "	return hash64shift( hash64shift( lo ) ^ hi );\n"; // perhaps something better LATER
		code += "}\n";
		code += "\n";
		code += "static float randof( float x, long long work_item, long long instance, long long step, int invocation_id ){\n";
		code += "	// Make a unique stamp for the random number sampled\n";
		code += "	// Unique factors: work item, tabular instance, serial number of RNG invocation in kernel, timestep \n";
		// TODO add a simulation properties digest, to decouple from exact timestep #
		// otherwise the time series for each individual random() invocation over time
		// 	would be an up/down scaled version of itself when fiddling with dt, instead of completely different
		// (do we want reproducibility or not, after all ??)
		code += "	// Capacities: 1T work items, 16M instances, 64K invocations, 1T timesteps \n";
		code += "	unsigned long long stamp_hi = work_item * (1ULL << 24) | instance % (1ULL << 24);\n";
		code += "	unsigned long long stamp_lo = invocation_id * (1ULL << 40) | step % (1ULL << 40);\n";
		code += "	unsigned long long sample = hash_128_to_64( stamp_hi, stamp_lo );\n";
		code += "	const/*ant*/int sample_scale = (1 << 23);\n";
		if(config.debug){
			code += "	printf(\"%llx\\n\", sample);\n";
		}
		code += "	float result = ( (float) ( sample % sample_scale ) ) / ( (float) (sample_scale) );\n";
		code += "	return x * result;\n";
		code += "}\n";
		code += "\n";
	};
	
	auto EmitWorkItemRoutineHeader = [ &config ]( std::string &code ){
		(void) config; // just in case
		code += "void doit( double time, float dt, const float *__restrict__ global_constants, long long const_local_index, \n"
		"const long long *__restrict__ global_const_table_f32_sizes, const Table_F32 *__restrict__ global_const_table_f32_arrays, long long table_cf32_local_index,\n"
		"const long long *__restrict__ global_const_table_i64_sizes, const Table_I64 *__restrict__ global_const_table_i64_arrays, long long table_ci64_local_index,\n"
		"const long long *__restrict__ global_state_table_f32_sizes, const Table_F32 *__restrict__ global_state_table_f32_arrays, Table_F32 *__restrict__ global_stateNext_table_f32_arrays, long long table_sf32_local_index,\n"
		"const long long *__restrict__ global_state_table_i64_sizes,       Table_I64 *__restrict__ global_state_table_i64_arrays, Table_I64 *__restrict__ global_stateNext_table_i64_arrays, long long table_si64_local_index,\n"
		"const float *__restrict__ global_state, float *__restrict__ global_stateNext, long long state_local_index, \n"
		"long long step ){\n";
		code +=   "	\n";
		
		
		code += "	\n";
		code += "	char initial_state = (step <= 0);\n";
		code += "	const float time_f32 = time; //when not accumulating small deltas, double precision is not necessary, and it messes up with SIMD\n";
		code +=   "	\n";
		
		code += "	const long long NOT_AN_INSTANCE = ~0xFee1600dLL; // if it's misused to index an array it will probably stop right there \xE3\x8B\xA1\n";
		code += "	long long instance = NOT_AN_INSTANCE; // for RNG use\n";
		code += "	long long rng_offset = 0; // for RNG use too\n";
		
		code += "	\n";
	};
	
	auto EmitWorkItemRoutineFooter = [ &config ]( std::string &code ){
		(void) config; // just in case
		
		code += "}\n";
	};
	auto EmitKernelFileFooter = [ &config ]( std::string &code ){
		(void) config; // just in case
		
		code += "// Generated code block END\n";
	};
	
	// LATER analyze cell types before generating codes, for compartment as work item
	printf("Implementing cell types...\n");
	// TODO build only the cells actually used
	for(size_t cell_seq = 0; cell_seq < cell_types.contents.size(); cell_seq++){
		const auto &cell_type = cell_types.contents[cell_seq];
		
		CellInternalSignature &sig = cell_sigs[cell_seq];
		
		sig.name = "Cell_type_"+itos(cell_seq);
		
		// TODO something more elegant to compile the actual code kernels:
		// dry run, or synchronization to use the same files(beware of distributed file systems !)
		#ifdef USE_MPI
		sig.name += "_rank_"+itos(my_mpi.rank);
		#endif
		
		printf("\nAnalyzing %s...:\n", sig.name.c_str());
		
		
		// cell-level work items for now
		SignatureAppender_Single AppendSingle_CellScope( sig.cell_wig );
		SignatureAppender_Table AppendMulti_CellScope( sig.cell_wig );
		InlineLems_AllocatorCoder DescribeLemsInline_CellScope( model, sig.cell_wig.random_call_counter, AppendSingle_CellScope, AppendMulti_CellScope );
		
		// standardize the nomenclature, yay!
		// <context>_<value or table>_<const, state, stateNext>
		// Exposed Context is a slice of the data space, that's relevant for a subitem (and which subsignatures refer to)
		// just like the work-item context is a slice of the global data space
		// Try not to over-use, to avoid overhead! Flatten wherever possible!
		// TODO make everything tabular, after all, and special-case flatten small vectors or sth
		auto ExposeSubitemContext = []( const std::string &to_context, const std::string &from_context, const std::string &tab ){
			std::string code;
			
			const auto &by = from_context, &to = to_context;
			
			code +=   "	const float *"+to+"_constants = "+by+"_constants + const_"+to+"_index;\n";
			code +=   "	const float *"+to+"_state     = "+by+"_state     + state_"+to+"_index;\n";
			code +=   "	      float *"+to+"_stateNext = "+by+"_stateNext + state_"+to+"_index;\n";
			code +=   "	\n";
				
			code += tab+"\tconst long long *"+to+"_const_table_f32_sizes      = "+by+"_const_table_f32_sizes      + table_cf32_"+to+"_index;\n";
			code += tab+"\tconst Table_F32 *"+to+"_const_table_f32_arrays     = "+by+"_const_table_f32_arrays     + table_cf32_"+to+"_index;\n";
			code += tab+"\tconst long long *"+to+"_const_table_i64_sizes      = "+by+"_const_table_i64_sizes      + table_ci64_"+to+"_index;\n";
			code += tab+"\tconst Table_I64 *"+to+"_const_table_i64_arrays     = "+by+"_const_table_i64_arrays     + table_ci64_"+to+"_index;\n";
			code += tab+"\tconst long long *"+to+"_state_table_f32_sizes      = "+by+"_state_table_f32_sizes      + table_sf32_"+to+"_index;\n";
			code += tab+"\tconst Table_F32 *"+to+"_state_table_f32_arrays     = "+by+"_state_table_f32_arrays     + table_sf32_"+to+"_index;\n";
			code += tab+"\t      Table_F32 *"+to+"_stateNext_table_f32_arrays = "+by+"_stateNext_table_f32_arrays + table_sf32_"+to+"_index;\n";
			code += tab+"\tconst long long *"+to+"_state_table_i64_sizes      = "+by+"_state_table_i64_sizes      + table_si64_"+to+"_index;\n";
			code += tab+"\t      Table_I64 *"+to+"_state_table_i64_arrays     = "+by+"_state_table_i64_arrays     + table_si64_"+to+"_index;\n";
			code += tab+"\t      Table_I64 *"+to+"_stateNext_table_i64_arrays = "+by+"_stateNext_table_i64_arrays + table_si64_"+to+"_index;\n";
				
			return code;
		};
		// utility function to access different contexts on same work item
		auto CloneSubitemIndices = []( const std::string &to_context, const std::string &from_context, const std::string &tab ){
			std::string code;
		
			const auto &by = from_context, &to = to_context;
			
			code +=   "	const long long const_"+to+"_index = const_"+by+"_index;\n";
			code +=   "	const long long state_"+to+"_index = state_"+by+"_index;\n";
			code +=   "	const long long table_cf32_"+to+"_index = table_cf32_"+by+"_index;\n";
			code +=   "	const long long table_ci64_"+to+"_index = table_ci64_"+by+"_index;\n";
			code +=   "	const long long table_sf32_"+to+"_index = table_sf32_"+by+"_index;\n";
			code +=   "	const long long table_si64_"+to+"_index = table_si64_"+by+"_index;\n";
			code +=   "	\n";
				
			return code;
		};
		
		
		// helpers for common cell parts (like synapses that may stick to any type of cell)
		
		// TODO use a crude approximation of Δi/Δv when conductance is missing!
		auto ExposeCurrentMaybeConductanceFromLems = [ ](const auto &comptype, const char *exponame_i, const char *exponame_g, const auto &tab){
			std::string ret;
			ret += tab+exponame_i+" = Lems_exposure_i;\n";
			if( comptype.common_exposures.conductance >= 0 ){
				ret += tab+exponame_g+" = Lems_exposure_g;\n";
			}
			return ret;
		};
				
		auto DescribeGenericSynapseInternals = [
			&model, &synaptic_components, &config, 
			ExposeCurrentMaybeConductanceFromLems
		](
			const std::string &tab, const std::string &for_what,
			const std::string &require_line, const std::string &expose_line,
			Int id_id,
			CellInternalSignature::SynapticComponentImplementation &synimpl,
			const SignatureAppender_Table &AppendMulti,
			const InlineLems_AllocatorCoder &DescribeLemsInline,
			std::string &internal_code
		){
			// use genericizable core implementations
			// for now, that means using the same expose/require line as LEMS
			// this is so synapse implementations (even blopla synapses) can be stuck in funny places like firing synapse inputs
			// if full optimization is desired, the genericizable core impl can be overriden for a synaptic-projection-specific implementation
			
			std::string &code = internal_code;
			char tmps[10000];
			
			const std::string Igap_suffix = Convert::Suffix(
				(Scales<Voltage>::native * Scales<Conductance>::native)
				.to(Scales<Current>::native)
			);
			const std::string Ichem_suffix = Igap_suffix;
			
			// XXX move these temp vars into the same closer scope where the exposure is done... or roll into ExposeCurrentMaybeConductanceFromLems
			code += tab+"	// Common exposures\n";
			code += tab+"	float Exposure_i = NAN;\n"; // the component *must* specify current
			code += tab+"	float Exposure_g = 0;\n"; // NOTE Exposing conductance is optional. If missing, assume 0 which reverts to Fwd Euler style handling 
			
			if(id_id < 0){
			
				SynapticComponent::Type core_id = SynapticComponent::Type(id_id + SynapticComponent::Type::MAX);
				
				code += tab+"	"+ require_line +"\n";
				
				switch(core_id){
				case SynapticComponent::Type::GAP :{
					const auto &for_that = for_what;
					std::string for_what = for_that + " Linear Gap Junction"; 
					// let's try derived constants
					code += "	// Linear gap junctions\n";
					
					size_t table_Gsyn  = synimpl.Table_Gbase = AppendMulti.Constant(for_what+" Base Conductance");
					
					sprintf(tmps, "		const float     *Gsyn_linear_gap  = local_const_table_f32_arrays[%zd];\n", table_Gsyn); code += tmps;
					
					sprintf(tmps, "		Exposure_i = Gsyn_linear_gap[instance] * (Vpeer - Vcomp)%s;\n", Igap_suffix.c_str()); code += tmps;
					// FIXME validate necessary conditions to aff Ggap to cable matrix
					sprintf(tmps, "		Exposure_g = Gsyn_linear_gap[instance] * 0;\n"); code += tmps;
					
					
					break;
				}
				case SynapticComponent::Type::EXP :{
					const auto &for_that = for_what;
					std::string for_what = for_that + " Exp Synapse";
					// gbase, erev, tau
					code += "	// Inbound exponential synapses\n";
					
					size_t table_Gbase = synimpl.Table_Gbase = AppendMulti.Constant(for_what+" Base Conductance");
					size_t table_Erev  = synimpl.Table_Erev  = AppendMulti.Constant(for_what+" Reversal Potential");
					size_t table_Tau   = synimpl.Table_Tau   = AppendMulti.Constant(for_what+" Time Constant");
					
					size_t table_G     = synimpl.Table_Grel  = AppendMulti.StateVariable(for_what+" Relative Conductance");
					
					sprintf(tmps, "	const float *Gbase_exp_one = local_const_table_f32_arrays[%zd];\n", table_Gbase); code += tmps;
					sprintf(tmps, "	const float *Erev_exp_one  = local_const_table_f32_arrays[%zd];\n", table_Erev); code += tmps;
					sprintf(tmps, "	const float *Tau_exp_one   = local_const_table_f32_arrays[%zd];\n", table_Tau); code += tmps;
					
					sprintf(tmps, "	const float *G_exp_one = local_state_table_f32_arrays[%zd];\n", table_G); code += tmps;
					
					sprintf(tmps, "	float   *Gnext_exp_one = local_stateNext_table_f32_arrays[%zd];\n", table_G); code += tmps;
					
					
					sprintf(tmps, "		Exposure_i = G_exp_one[instance] * ( Erev_exp_one[instance] - Vcomp)%s;\n", Ichem_suffix.c_str());  code += tmps;
					sprintf(tmps, "		Exposure_g = G_exp_one[instance] ;\n");  code += tmps;
					
					code   += "		if(!initial_state){\n";
					sprintf(tmps, "			Gnext_exp_one[instance] = G_exp_one[instance] - dt * ( G_exp_one[instance] / Tau_exp_one[instance] )%s;\n", "");  code += tmps;
					code   += "		}else{\n";
					sprintf(tmps, "			Gnext_exp_one[instance] = G_exp_one[instance];");  code += tmps;
					code   += "		}\n";
					// NOTE synapse state has been initialized from synapse default value, otherwise it would have been initialized here
					
					// TODO keep these loops separated, to support lazy synapse triggering
					code   += tab+"	if(!initial_state){\n";
					code   += tab+"		if( spike_in_flag ) {\n";
					if(config.debug){
					code   += tab+"			printf(\"kaboom, baby! %lld\\n\", instance);\n";
					}
					code   += tab+"			Gnext_exp_one[instance] = G_exp_one[instance] + Gbase_exp_one[instance];\n";
					code   += tab+"		}\n";
					code   += tab+"	}\n";
					
					break;
				}
			
				default:
					// internal error
					printf("internal error: Unknown synaptic component core_id %d\n", core_id);
					return false;
				}
				
				code += tab+expose_line+"\n";
			}
			else{
				// Not just LEMS, but any custom type
				// id_id cannot aggregate structurally different components of same "type" (such as blocking plastic )
				// TODO abolish id_id and use simple case impl'd flag instead
				// for core types, and identical component (conflicts with constant deduplication TODO)
				const auto &for_that = for_what;
				// TODO more details in this description 
				
				Int syncomp_seq = id_id; // TODO component type
				const auto &syncomp = synaptic_components.get(syncomp_seq);
				// TODO eliminate redundant allocations for same component type
				
				code += "	{\n"; // special handling start
				auto ExposeFromLems = [&](const auto &comptype){ return ExposeCurrentMaybeConductanceFromLems(comptype, "Exposure_i", "Exposure_g", tab+"\t"); };
				// now what to do with the special type
				if(syncomp.type == SynapticComponent::Type::BLOCKING_PLASTIC ){
					std::string for_what = for_that + " Blocking/Plastic Synapse";
					
					const auto &blopla_inst = syncomp.component;
					const auto &blopla_type = model.component_types.get(blopla_inst.id_seq);
					auto &blopla_subsig = synimpl.synapse_component;
					
					
					blopla_subsig = DescribeLems::AllocateSignature(blopla_type, blopla_inst, &AppendMulti, for_what + " Component LEMS");
					
					
					// fuse dynamics of three LEMS elements, as below:
					//code += DescribeLemsInline.TableLoop( tab, for_what, blopla_subsig );
					//code += tab+"{\n";
						
						// expose requirements, so they can be shared with child sub-components
						code += tab+require_line+"\n";
						// first compute the sub-components to expose block, plasticity factors
						code += tab+"	float block_factor = 1, plasticity_factor = 1;";
						
						// NOTE assume blocking/plastic are lemsified, if they exist
						if( syncomp.blopla.block_mechanism.type != SynapticComponent::BlockingPlasticSynapse::BlockMechanism::NONE ){
							const auto &blo_inst = syncomp.blopla.block_mechanism.component;
							const auto &blo_type = model.component_types.get(blo_inst.id_seq);
							auto &blo_subsig = synimpl.block_component;
							blo_subsig = DescribeLems::AllocateSignature(blo_type, blo_inst, &AppendMulti, for_what + " Block Component");
							code += tab+"{\n";
							code += DescribeLemsInline.TableInner( tab, for_what + " Block Component", blo_type, blo_subsig, "", "block_factor = Lems_exposure_blockFactor;", config.debug );
							code += tab+"}\n";
						}
						
						if( syncomp.blopla.plasticity_mechanism.type != SynapticComponent::BlockingPlasticSynapse::PlasticityMechanism::NONE ){
							const auto &pla_inst = syncomp.blopla.plasticity_mechanism.component;
							const auto &pla_type = model.component_types.get(pla_inst.id_seq);
							auto &pla_subsig = synimpl.plasticity_component;
							pla_subsig = DescribeLems::AllocateSignature(pla_type, pla_inst, &AppendMulti, for_what + " Plasticity Component");
							code += tab+"{\n";
							code += DescribeLemsInline.TableInner( tab, for_what + " Plasticity Component", pla_type, pla_subsig, "", "plasticity_factor = Lems_exposure_plasticityFactor;", config.debug );
							code += tab+"}\n";
						}
						
						// then compute the synaptic component
						code += DescribeLemsInline.TableInner( tab, for_what, blopla_type, blopla_subsig, "", ExposeFromLems(blopla_type)+expose_line, config.debug );
						
					//code += tab+"}\n";
					
				}
				else if( syncomp.component.ok() ){
					std::string for_what = for_that + " LEMS Synaptic Component";
					
					// allow loop to be overridden externally; don't use the all-in-one solution that includes its own loop
					const ComponentInstance &compinst = syncomp.component;
					const auto &comptype = model.component_types.get(compinst.id_seq);
					CellInternalSignature::ComponentSubSignature &compsubsig = synimpl.synapse_component;
					
					synimpl.synapse_component = DescribeLems::AllocateSignature(comptype, compinst, &AppendMulti, for_what);
					code += DescribeLemsInline.TableInner( tab+"\t", for_what, comptype, compsubsig, require_line, ExposeFromLems(comptype)+expose_line, config.debug );
				}
				else{
					printf("internal error: synaptic component %ld is neither special case nor lemsified \n", syncomp_seq);
					return false;
				}
				
				
				code += tab+"	}\n"; // special handling end
			}
			
			return true;
		};
		
		// generate tables and code for each synapse type
		auto ImplementSynapseType = [ &DescribeGenericSynapseInternals, &model, &config, &synaptic_components ](
			const SignatureAppender_Single &AppendSingle, const SignatureAppender_Table &AppendMulti,
			const InlineLems_AllocatorCoder &DescribeLemsInline,
			Int &random_call_counter,
			const std::string &for_what,
			const std::string &tab,
			Int id_id,
			auto &synapse_impls,
			std::string &ccde
		){
			char tmps[1000];
			ccde += tab+"{\n";
			
			CellInternalSignature::SynapticComponentImplementation synimpl;
			
			// TODO this is a bit dirty, refactor by relegating id_id usage to final deduplication, or sth
			SynapticComponent fake_syn;
			if( id_id < 0 ){
				fake_syn.type = SynapticComponent::Type(id_id + SynapticComponent::Type::MAX);
			}
			else fake_syn = synaptic_components.get(id_id);
			
			const bool needs_spike = fake_syn.HasSpikeIn(model.component_types);
			const bool needs_Vpeer = fake_syn.HasVpeer(model.component_types);
		
			// add weight for all synapse types, on the same side as the realized syn. component (ie post for chemical)
			size_t table_Weight = synimpl.Table_Weight  = AppendMulti.Constant(for_what+" Weight");

			sprintf(tmps, "	const float     *Weight  = local_const_table_f32_arrays[%zd];\n", table_Weight); ccde += tmps;
			
			// allocate and expose tables for spike and Vpeer communication, as needed
			if( needs_spike ){
				
				size_t table_Trig     = synimpl.Table_Trig  = AppendMulti.StateI64(for_what+" Trigger");
				sprintf(tmps, "	long long   *Trigger = local_state_table_i64_arrays[%zd];\n", table_Trig); ccde += tmps;
				
				bool uses_delay = true; // maybe elide LATER
				if( uses_delay ){
					size_t Table_Delay = synimpl.Table_Delay  = AppendMulti.Constant(for_what+" Delay");
					sprintf(tmps, "	const float *Delay = local_const_table_f32_arrays[%zd];\n", Table_Delay); ccde += tmps;
					
					size_t Table_NextSpike = synimpl.Table_NextSpike  = AppendMulti.StateVariable(for_what+" Next Spike");
					sprintf(tmps, "	const float *NextSpike = local_state_table_f32_arrays[%zd];\n", Table_NextSpike); ccde += tmps;
					sprintf(tmps, "	float *NextSpike_Next = local_stateNext_table_f32_arrays[%zd];\n", Table_NextSpike); ccde += tmps;
				}
				
			}
			if( needs_Vpeer ){
				// LATER The tricky part is to get both GJ components to use the same weight array LATER? probably not, due to access irregularities?
				
				size_t table_Vpeer = synimpl.Table_Vpeer = AppendMulti.ConstI64(for_what+" Vpeer Global State Index");
				
				sprintf(tmps, "	const long long *Vpeer_array = local_const_table_i64_arrays[%zd];\n", table_Vpeer); ccde += tmps;
				
			}
			
			// and common parts for all syn.components
			
			// expose instances
			if( needs_spike ){
				sprintf(tmps, "	const long long Instances = local_state_table_i64_sizes[%zd]; //same for all parallel arrays\n", synimpl.Table_Trig ); ccde += tmps;
			}
			else if( needs_Vpeer ){
				sprintf(tmps, "	const long long Instances = local_const_table_i64_sizes[%zd]; //same for all parallel arrays\n", synimpl.Table_Vpeer); ccde += tmps;
			}
			else if( id_id >= 0 && fake_syn.component.ok() ){
				ccde += DescribeLemsInline.TableInstances( tab, for_what + "LEMS Component", synimpl.synapse_component );
			}
			else{
				// TODO something LEMS-like? for generality?
				// LEMS handles it itself, so put it here
				// otherwise, another way shall be found LATER
				printf("internal error: synapse type %ld should receive spikes or Vpeer, or have LEMS properties, or any other way to determine its physical existence\n", id_id);
				return false;
			}
			
			auto GetGenericRequirementCode = [ & ](){
				std::string require_line;
				
				bool uses_weight = true; // maybe elide LATER
				if( uses_weight ){
					require_line += "\n" + tab + "float weight = Weight[instance];\n";
				}
				
				
				if( needs_spike ){
					bool uses_delay = true; // maybe elide LATER
					
					require_line += "\n" + tab +
						// NB: this part should be modified along with pre-synaptic spike sending !
						"char spike_in_flag = 0;\n"
						+tab+"if( !initial_state ){\n" // TODO mark initial_state as rarely taken
						+tab+"\tspike_in_flag = !!Trigger[instance];\n" // TODO mask
						+tab+"\tTrigger[instance] = 0;\n" // XXX refactor clearing the flag somewhere else, to avoid accidentally reading it as zero in multiple inclusions!
						+tab+"}\n"
					; // LATER use mask
					
					if( uses_delay ){
						require_line += "float delay = Delay[instance];\n"
						+tab+ "float next_spike = NextSpike[instance];\n"
						+tab+ "float next_next_spike = next_spike;\n"
						+tab+ "char spike_now = 0;\n"
						+tab+"if( !initial_state ){\n" // TODO mark initial_state as rarely taken
						
						+tab+"if( time_f32 <= next_spike && next_spike < time_f32 + dt ){\n"
						+tab+"	spike_now = 1;\n"
						+tab+"}\n"
						
						+tab+"if( spike_in_flag ){\n"
						+tab+"	float fresh_spike = time_f32 + delay;\n"
						+tab+"	if( time_f32 <= fresh_spike && fresh_spike < time_f32 + dt ){\n"
						+tab+"		spike_now = 1;\n"
						+tab+"	}\n"
						+tab+"	if( next_next_spike < time_f32 + dt && next_next_spike < fresh_spike ){\n"
						+tab+"		next_next_spike = fresh_spike; // keep first incoming spike\n"
						+tab+"	}\n"						
						+tab+"}\n"
						
						+tab+"}else{\n"
						
						+tab+"}\n"
						
						+tab+"spike_in_flag = spike_now;\n"
						+tab+"NextSpike_Next[instance] = next_next_spike;\n"
						;
					}
				}
				
				if( needs_Vpeer ){
					// TODO refactor better, possibly re-use the global_state vector to improve performance
					// require_line += "\n" + tab + "float Vpeer = global_state[ Vpeer_array[instance] ];";
					// TODO make this an inline function
					require_line += "\n" + tab + "float Vpeer;";
					require_line += "\n" + tab + "{";
					require_line += "\n" + tab + "\tconst unsigned long long packed_id =  Vpeer_array[instance];";
					require_line += "\n" + tab + "\tconst unsigned long long table_id = packed_id / (1 << 24);";
					require_line += "\n" + tab + "\tconst unsigned long long entry_id = packed_id % (1 << 24);";
					if(config.debug){
					require_line += "\n" + tab + "\tprintf(\"vpe %llx\\t%llu\\t%llu\\t%p\\n\", packed_id, table_id, entry_id,global_state_table_f32_arrays[table_id]);";
					require_line += "\n" + tab + "\tfflush(stdout);";
					}
					require_line += "\n" + tab + "\tVpeer = global_state_table_f32_arrays[table_id][entry_id];";
					// require_line += "\n" + tab + "\tprintf(\"vpeo.\\n\");";
					// require_line += "\n" + tab + "\tfflush(stdout);";
					
					require_line += "\n" + tab + "}";
				}
				
				return require_line;
			};
			auto GetGenericExposureCode = [ & ](){
				
				return std::string("I_syn_aggregate += Exposure_i * weight; G_syn_aggregate += Exposure_g * weight;\n"); // TODO merge this line with input code as well
			};
			
			const std::string expose_line = GetGenericExposureCode();
			
			std::string require_line = GetGenericRequirementCode();
			
			
			ccde += tab+"float I_syn_aggregate = 0;\n";
			ccde += tab+"float G_syn_aggregate = 0;\n";
			
			// Simple not-so-scalable solution:
			// scan everything, don't use lazy triggering
			ccde   += tab+"for(long long instance = 0; instance < Instances; instance++){\n";
			
			
			
			std::string syn_internal_code;
			if( !DescribeGenericSynapseInternals( tab, for_what, require_line, expose_line, id_id, synimpl, AppendMulti, DescribeLemsInline, syn_internal_code ) ) return false;
			ccde   += syn_internal_code;
			
			ccde   += tab+"}\n"; // for loop end
			// add the gathered curent to total synapse current
			sprintf(tmps, "	I_synapses_total += I_syn_aggregate;\n"); ccde += tmps;
			sprintf(tmps, "	G_synapses_total += G_syn_aggregate;\n"); ccde += tmps;
			ccde   += "\n";
			
			synapse_impls[id_id] = synimpl; // gap junctions and chemical synapses have different id_id's
			
			ccde += tab+"}\n";
			
			return true;
		};
		auto ImplementInputSource = [ &DescribeGenericSynapseInternals, &GetSynapseIdId, &ExposeCurrentMaybeConductanceFromLems, &config, &model, &input_sources ](
			const SignatureAppender_Single &AppendSingle, const SignatureAppender_Table &AppendMulti,
			const InlineLems_AllocatorCoder &DescribeLemsInline,
			Int &random_call_counter,
			const std::string &for_what,
			const std::string &tab,
			Int id_id,
			auto &input_impls,
			std::string &ccde
		){
			char tmps[1000];
			CellInternalSignature::InputImplementation inpimpl;
			ccde += tab+"{\n";
			
			inpimpl.Table_Weight     = AppendMulti.Constant(for_what+" Weight");
			sprintf(tmps, "	const float     *Weight     = local_const_table_f32_arrays[%zd];\n", inpimpl.Table_Weight); ccde += tmps;
			
			// helpers
			auto ImplementTabularSpikeList_OpenEnd = [ &config, &AppendMulti, &tmps ]( const std::string &for_what, const std::string &tab, auto &inpimpl, auto &ccde){
				(void) config; // just in case
				
				// for each instance of the same input source type:
				// a slice of a common spike time vector, and a start index to begin from for each instance
				// and a state variable to which index is coming up for each instance
				size_t table_Times = inpimpl.Table_SpikeListTimes  = AppendMulti.Constant(for_what+" Spike Times");
				//size_t table_Start = inpimpl.Table_SpikeListStarts = AppendMulti.ConstI64(for_what+" Spike Index Start");
				size_t table_Posit = inpimpl.Table_SpikeListPos    = AppendMulti.StateI64(for_what+" Spike Index Position");
				
				// positions are initialized at initial tables time, yay!
	
				sprintf(tmps, "	const long long Instances = local_state_table_i64_sizes[%zd]; //same for all parallel arrays\n", inpimpl.Table_SpikeListPos ); ccde += tab+tmps;
				
				ccde   += tab+"for(long long instance = 0; instance < Instances; instance++){\n";
			
				sprintf(tmps, "const float     *Spike_Times  = local_const_table_f32_arrays[%zd];\n", table_Times); ccde += tab+tmps;
				sprintf(tmps, "const long long *Positions  = local_state_table_i64_arrays[%zd];\n", table_Posit); ccde += tab+tmps;
				sprintf(tmps, "      long long *PositNext  = local_stateNext_table_i64_arrays[%zd];\n", table_Posit); ccde += tab+tmps;
				
				// TODO wrap into a reqstring ?
				ccde   += tab+"char spiker_fired_flag = 0;\n";
				ccde   += tab+"long long pos = Positions[instance];\n";
				ccde   += tab+"while( time_f32 >= Spike_Times[pos] ){\n";
				ccde   += tab+"	spiker_fired_flag = 1;\n";
				ccde   += tab+"	pos++;\n";
				ccde   += tab+"}\n";
				ccde   += tab+"if( !initial_state ){\n";
				ccde   += tab+"	PositNext[instance] = pos;\n";
				ccde   += tab+"}\n";
				
				return true;
			};
			
			if(id_id < 0){
				InputSource::Type core_id = InputSource::Type(id_id + InputSource::Type::MAX);
				
				switch(core_id){
					
				case InputSource::Type::PULSE :{
					const auto &for_that = for_what; // because C won't recognize the outer scope variable in initialization >:(
					std::string for_what = for_that + " DC Pulse";
					// let's try derived constants
					ccde += "	// Pulse inputs\n";
					
					size_t table_Imax     = inpimpl.Table_Imax     = AppendMulti.Constant(for_what+" Imax");
					size_t table_start    = inpimpl.Table_Delay    = AppendMulti.Constant(for_what+" Start");
					size_t table_duration = inpimpl.Table_Duration = AppendMulti.Constant(for_what+" Duration");
					
					sprintf(tmps, "	const long long Instances_input_pulse = local_const_table_f32_sizes[%zd]; //same for all parallel arrays\n", table_Imax); ccde += tmps;
					
					sprintf(tmps, "	const float     *Imax_input_pulse     = local_const_table_f32_arrays[%zd];\n", table_Imax); ccde += tmps;
					sprintf(tmps, "	const float     *Start_input_pulse    = local_const_table_f32_arrays[%zd];\n", table_start); ccde += tmps;
					sprintf(tmps, "	const float     *Duration_input_pulse = local_const_table_f32_arrays[%zd];\n", table_duration); ccde += tmps;
					
					sprintf(tmps, "	float I_input_pulse = 0;\n"); ccde += tmps;
					if(config.use_icc){
						ccde +=   "	 #pragma novector\n";
					}
					sprintf(tmps, "	for(long long instance = 0; instance < Instances_input_pulse; instance++){\n"); ccde += tmps;
					ccde +=   "		if( Start_input_pulse[instance] <= time && time <=  Start_input_pulse[instance] +  Duration_input_pulse[instance] ) I_input_pulse += Imax_input_pulse[instance] * Weight[instance];\n";
					ccde   += "	}\n";
					sprintf(tmps, "	I_input_total += I_input_pulse;\n"); ccde += tmps;
					sprintf(tmps, "	I_input_total += 0;\n"); ccde += tmps; // the ideal current probe is invariant with voltage
					ccde   += "\n";
					
					break;
				}
				case InputSource::Type::SPIKE_LIST :{
					const auto &for_that = for_what; // because C won't recognize the outer scope variable in initialization >:(
					std::string for_what = for_that + " Spike List";
					
					if( !ImplementTabularSpikeList_OpenEnd( for_what, tab, inpimpl, ccde ) ) return false;
					
					ccde   += tab+"spike_in_flag |= spiker_fired_flag;\n";
					ccde   += tab+"}\n"; // for loop end. TODO more elegant
					
					break;
				}
				default:
					// internal error
					printf("Unknown input core_id %d\n", core_id);
					return false;
				}
			}
			else{
				// Not just LEMS, but any custom type
				const auto &for_that = for_what;
				
				Int input_source_seq = id_id;
				const auto &input_source = input_sources.get(id_id);
				// TODO eliminate redundant allocations for same component type
				
				std::string require_line = "float weight = Weight[instance];";
				
				ccde += tab+"{\n";
				ccde += tab+"float I_syn_aggregate = 0;\n";
				ccde += tab+"float G_syn_aggregate = 0;\n";
				
				
				if(
					input_source.type == InputSource::Type::TIMED_SYNAPTIC
					|| input_source.type == InputSource::Type::POISSON_SYNAPSE
					|| input_source.type == InputSource::Type::POISSON_SYNAPSE_TRANSIENT
				){
					std::string for_what;
					
					// de-dupe by (synapse type, coretype) LATER
					// first implement the spiker
					if( input_source.type == InputSource::Type::TIMED_SYNAPTIC ){
						for_what = for_that + " Timed Synaptic Input";
						
						if( !ImplementTabularSpikeList_OpenEnd( for_what, tab, inpimpl, ccde ) ) return false;
					}
					else if(
						input_source.type == InputSource::Type::POISSON_SYNAPSE
						|| input_source.type == InputSource::Type::POISSON_SYNAPSE_TRANSIENT
					){
						for_what = for_that + " Poisson Firing Synapse";
						
						const auto &spik_inst = input_source.component;
						const auto &spik_type = model.component_types.get(spik_inst.id_seq);
						auto &spik_subsig = inpimpl.component;
						spik_subsig = DescribeLems::AllocateSignature(spik_type, spik_inst, &AppendMulti, for_what + " Spiker");
						
						ccde += DescribeLemsInline.TableLoop( tab, for_what, spik_subsig );
						ccde += tab+"{\n";
						
						ccde   += tab+"char spiker_fired_flag = 0;\n";
						
						ccde   += tab+"{\n"; // spiker start
						ccde += DescribeLemsInline.TableInner( tab, for_what + " Spiker", spik_type, spik_subsig, "", "spiker_fired_flag = Lems_eventout_spike;", config.debug );
						ccde   += tab+"}\n"; // spiker end
					}
					else{
						printf("internal error: input component %ld code for what sort of firing synapse input? \n", input_source_seq);
						return false;
					}
					
					ccde   += tab+"char spike_in_flag = spiker_fired_flag;\n";
					
					// weight is applied here since each synapse is an input instance at the same time
					std::string expose_line = "I_syn_aggregate += Exposure_i * weight; G_syn_aggregate += Exposure_g * weight;";
					
					std::string syn_internal_code;
					if( !DescribeGenericSynapseInternals( tab, for_what, require_line, expose_line, GetSynapseIdId(input_source.synapse), inpimpl.synimpl, AppendMulti, DescribeLemsInline, syn_internal_code ) ) return false;
					ccde   += syn_internal_code;
					
					ccde   += tab+"}\n"; // for loop end
					
					sprintf(tmps, "I_input_total += I_syn_aggregate;\n"); ccde += tab+tmps;
					sprintf(tmps, "G_input_total += G_syn_aggregate;\n"); ccde += tab+tmps;
					ccde   += tab+"\n";
					
				}
				else if( input_source.component.ok() ){
	
					std::string for_what = for_that + " LEMS Input";
					
					const auto &comptype = model.component_types.get(input_source.component.id_seq);
				
					// use temp variables in the place of lems exposures because g might not be exposed
					// It's done at the scope of component based inputs because hardcoded inputs handle the loop in their inner code, unlike synpases which are fully decoupled.
					// TODO decouple input internals for composite inputs
					// TODO aggregate seprately for each input type, instead of accumulating on input total
					std::string expose_line;
					// XXX move these temp vars into the same closer scope where the exposure is done... or roll into ExposeCurrentMaybeConductanceFromLems
					expose_line += tab+"	float Exposure_i = NAN;\n"; // the component *must* specify current
					expose_line += tab+"	float Exposure_g = 0;\n"; // NOTE Exposing conductance is optional. If missing, assume 0 which reverts to Fwd Euler style handling 
					expose_line += ExposeCurrentMaybeConductanceFromLems(comptype, "Exposure_i", "Exposure_g", tab);
					expose_line += "	I_input_total += Exposure_i * weight; G_input_total += Exposure_g * weight;\n";
				
					// TODO return bool for error handling
					ccde += DescribeLemsInline.TableFull( input_source.component, "\t", for_what, inpimpl.component, require_line, expose_line, config.debug );
				}
				else{
					printf("internal error: input component %ld is neither special case nor lemsified \n", input_source_seq);
					return false;
				}
				
				ccde += "	}\n"; // special case end
			}
			
			ccde += tab+"}\n";
			input_impls[id_id] = inpimpl;
			return true;
		};
		auto ImplementSpikeSender = [ &config ](
			const std::string &condition,
			const SignatureAppender_Table &AppendMulti,
			const std::string &for_what,
			CellInternalSignature::SpikeSendingImplementation &spiker,
			std::string &code
		){
			char tmps[1000];
			// one table, with  destination indexes
			// pack 1 trillion tables -> 16 million entries into 64bit indexes, upgrade if needed LATER
			
			size_t table_Spike_recipients = spiker.Table_SpikeRecipients = AppendMulti.ConstI64( for_what+" Spike recipients");
			// printf("spiker send %zd %zd\n", cell_seq, comp_seq);
			
			sprintf(tmps, "	const long long Instances_Spike_recipients = local_const_table_i64_sizes[%zd]; //same for all parallel arrays\n", table_Spike_recipients); code += tmps;
			sprintf(tmps, "	const long long *Spike_recipients          = local_const_table_i64_arrays[%zd];\n", table_Spike_recipients); code += tmps;
			// NOTE Vnext should have been calculated by this point!
			code   += "	// Spike check\n";
			code   += "	if( " + condition + " ) {\n";
			code   += "		for(long long instance = 0; instance < Instances_Spike_recipients; instance++){\n";
			code   += "			const unsigned long long packed_id = Spike_recipients[instance];\n";
			code   += "			const unsigned long long table_id = packed_id / (1 << 24);\n";
			code   += "			const unsigned long long entry_id = packed_id % (1 << 24);\n";
			code   += "			const unsigned long long word_id = entry_id >> 0;\n";
			code   += "			const unsigned long long mask = 1 << 0;\n";
			if(config.debug){
				code   += "			printf(\"%p %p %llx %llu %llu %llu\\n\", global_stateNext_table_i64_arrays, global_stateNext_table_i64_arrays[table_id], packed_id, table_id, entry_id, word_id);\n";
			}
			
			
			//TODO what happens if a spike occurs while the synapse is already active ?? reset? accumulate? ( what sort of accumulation ? )
			// NeuroML does not offer a way for two different spike sources to trigger the same post-synaptic instances,
			// which means consecutive spikes on the same post-synaptic instance originate from the same spike source instance, with probably the same delay.
			// If multiple spikes must concurrently travel along a pure lag component, with each spike remaining a distinct time of occurrence
			// (such as the case of new spikes accumulating on an already active synapse), possible solutions are:
			// 	- discretize the medium of propagation through multiple compartments 
			// 	- implement a general spike priority queue, or specialized components
			// 		- parallel priority queues do exist, even for gpu's; may need some modifications for up-to-current-time popping, though (or just reinsert what is not yet ready)
			// 	- assume a finite number of concurrently propagating spikes (and fall back to general priority queue otherwise)
			
			//sprintf(tmps, "			global_stateNext_table_i64_arrays[table_id][word_id] = mask;\n" );  code += tmps;
			//sprintf(tmps, "			atomic_fetch_or_explicit( (atomic_ullong *) &( global_stateNext_table_i64_arrays[table_id][word_id] ), mask, memory_order_relaxed );\n" );  code += tmps;
			code   += "			__sync_fetch_and_or( &( global_stateNext_table_i64_arrays[table_id][word_id] ), mask );\n" ;
			
			// end spike sending case
			code   += "		}\n";
			code   += "	}\n";
			
			return true;
		};
		
		auto ImplementRngSeed = [ &config, &model ](
			const SignatureAppender_Single &AppendSingle,
			const std::string &for_what, const std::string &tab,
			const std::string &subitem_context,
			auto &rng_impl,
			std::string &ccde
		){
			char tmps[1000];
			CellInternalSignature::InputImplementation inpimpl;
			
			ptrdiff_t Index_RngSeed = rng_impl.Index_RngSeed = AppendSingle.Constant( 0, for_what+" Cell RNG Seed" );
			
			sprintf(tmps, "const int cell_rng_seed = EncodeF32ToI32(%s_constants[%zd]);\n", subitem_context.c_str(), Index_RngSeed); ccde += tab+tmps;
			
			return true;
		};
		
		// generate per-cell data signature, initial data and iteration code
		
		if( cell_type.type == CellType::PHYSICAL ){
		
		const PhysicalCell &cell = cell_type.physical;
		const Morphology &morph = morphologies.get(cell.morphology);
		const BiophysicalProperties &bioph = biophysics.at(cell.biophysicalProperties);
		
		CellInternalSignature::PhysicalCell &pig = sig.physical_cell;
		const auto &comp_disc = pig.compartment_discretization;
		const int32_t number_of_compartments = (int32_t) comp_disc.number_of_compartments;
		
		const ScaleEntry microns = {"um", -6, 1.0}; // In NeuroML, Morphology is given in microns
		
		// preprocess Morphology for geometry, connectivity info, and compartmental subdivision
		// interaction between state variables is, of course, also determined by connectivity between compartments
		// NeuroML assumes a tree model of truncated cone-shaped 'segments'
		// 	( which may be further subdivided in compartments, for compatibility with NEURON )
		std::vector< std::vector<int32_t> > comp_connections(number_of_compartments);
		
		// physical properties
		std::vector<float> compartment_lengths(number_of_compartments, NAN);
		std::vector<float> compartment_areas  (number_of_compartments, NAN);
		std::vector<float> compartment_volumes(number_of_compartments, NAN);
		std::vector<float> compartment_path_length_from_root(number_of_compartments, NAN);
		
		
		// TODO:
		// If the values of properties have to be sampled across the neuron, there are two ways to go about it:
		// The first one is to sample only for the middles of compartments. 
		// This is more accurate if few segments are split into many compartments, and less accurate if many segments are merged into few compartments.
		// It has the additional benefit that all calculations are performed compartment-wise, which is amenable to in-line initialisation in the per-compartment code.
		// The second one is to sample for the middles of segments.
		// This is less accurate if few segments are split into many compartments, and more accurate if many segments are merged into few compartments.
		// It has the complication of mapping segment values to compartment values, through the conversion process of said segments into compartments; which is not as amenable to parallel run-time code.
		// A third one is to apply all segment-and compartment-wise intersections, and evaluate each "compartment fragment" that comes out of that.
		// All these approaches assume that said values can be interpolated smoothly along the length of neurites.
		
		// Another problem comes up with values which are simply not integrable over the extent of a neurite.
		//  Examples are starting membrane voltage, and spike initiation "threshold" voltage. Exclusively one value can hold for these oer a compartment.
		
		// For now:
		// - The geometric properties of the neurite compartments (area, volume) are precisely integrated using the linearly varying per segment values. This is convenient because these geometric values are laid out explicitly.
		// - The resistance of the neurite is also precisely evaluated through the per segment values. To do otherwise would be to calculate the (1/length) factor to multiply by the per-compartment resistivity.
		// - The conductance of density mechanisms is *IM*precisely evaluated by mapping the per-segment parameters over the compartments which overlap with the segments.
		// - The conductance of inhomogeneous ion channel distributions is *IM*precisely evaluated through the per compartment central path lengths. This is possible through the pre-calculated per-compartment area.
		
		
		
		// Density mechanisms such as ion channels and ion pools are applied over the whole compartments that the designated segment group to apply over overlaps with, whether the overlap is partial or not.
		// A possible unintended consequence is for a compartment that is straddling two segments with different mechanism specifications over the two segments, to contain both sets of mechanisms, each with excess conductance.
		// This could be refined by introducing another paameter of "effective area over a compartment" to the numerical model of each mechanicm instance. 
		// A hacky alternative could be to re-scale the Gmax, or ρmax, but this only applies if the effect is a linear function of area. (Also does not necessarily work with LEMS mechs.)
		// But in principle, the proper way to ensure that mechnaism are applied to a specific compartment is to make an explicit semgment to represent the compartment.
		
		// process connectivity
		printf("\tAnalyzing internal connectivity...\n");
		for( Int comp_seq = 0; comp_seq < number_of_compartments; comp_seq++ ){
			const auto comp_parent = comp_disc.tree_parent_per_compartment[comp_seq];
			
			if(!(comp_parent < 0)){
				comp_connections[comp_seq].push_back(comp_parent);
				comp_connections[comp_parent].push_back(comp_seq);
			}
			// Since comp_parent < comp_seq, each comp_connections[i] list is always going to be ordered !
		}
		
		// process geometry of morphology
		printf("\tAnalyzing geometry...\n");
		for( int comp_seq = 0; comp_seq < number_of_compartments; comp_seq++ ){
			
			auto &comp_length = compartment_lengths[comp_seq];
			auto &comp_area   = compartment_areas  [comp_seq];
			auto &comp_volume = compartment_volumes[comp_seq];
			
			// summate by iterating through all segment-based fractions
			comp_length = 0;
			comp_area   = 0;
			comp_volume = 0;
			
			const auto &segs_in_comp = comp_disc.compartment_to_segment_seq[comp_seq];
			
			for( int seg_in_comp = 0; seg_in_comp < (int) segs_in_comp.size(); seg_in_comp++ ){
				Real from_fractionAlong = (seg_in_comp == 0) ? comp_disc.compartment_to_first_segment_start_fractionAlong[comp_seq] : 0;
				Real   to_fractionAlong = (seg_in_comp+1 == (int) segs_in_comp.size()) ? comp_disc.compartment_to_last_segment_end_fractionAlong[comp_seq] : 1;
				Int seg_seq = segs_in_comp[seg_in_comp];
				
				// TODO check if there is any case proximal should be the specified one, or the distal of parent seg instead
				
				const Morphology::Segment &seg = morph.segments.atSeq(seg_seq);
				
				// NOTE the actual 3D points used may differ from what's stated in the NeuroML tag
				// following obscure rules followed by the NeuroML exporter
				// TODO deduplicate the intra-seg mapping through lerp, to make sure it's consistent
				
				const Morphology::Segment::Point3DWithDiam seg_from = LerpPt3d(seg.proximal, seg.distal, from_fractionAlong);
				const Morphology::Segment::Point3DWithDiam seg_to   = LerpPt3d(seg.proximal, seg.distal,   to_fractionAlong);
				
				Real comp_fragment_length = DistancePt3d( seg_from, seg_to );
				
				comp_length += comp_fragment_length;
				comp_area   += GeomHelp::Area( comp_fragment_length, seg_from.d, seg_to.d );
				comp_volume += GeomHelp::Volume( comp_fragment_length, seg_from.d, seg_to.d );
				
			}
		}
		// get path length from root too
		// comp_seq order puts parents first always, just like seg_seq order
		for( int comp_seq = 0; comp_seq < number_of_compartments; comp_seq++ ){
			const auto comp_parent = comp_disc.tree_parent_per_compartment[comp_seq];
			
			if( comp_parent < 0 ){
					// well, Neuron assumes the *proximal edge* of a compartment to be the beginning of "Path length from root",
					// (thus diastnnce(comp0) = comp0_length/2) rather than the middle of the compartment (in which case it would be disctance(comp0) = 0).
					// compartment_path_length_from_root[comp_seq] = 0;
					compartment_path_length_from_root[comp_seq] = (compartment_lengths[comp_seq] / 2);
					
			}
			else{
				compartment_path_length_from_root[comp_seq] = 
					compartment_path_length_from_root[comp_parent] 
					+ (compartment_lengths[comp_parent] + compartment_lengths[comp_seq])/2 ;
			}
		}
		
		// process passive biophysics
		printf("\tAnalyzing cable equation...\n");
		// membrane specific capacitance, axial resistivity for every segment, sample them to integrate capacitance and resistance for each compartment.		
		// d_lambda rule can also be calculated from Cm and Ra (NEURON book, chapter 5)
		std::vector<Real> segment_Cm_(morph.segments.contents.size(), NAN);
		std::vector<Real> segment_Ra_(morph.segments.contents.size(), NAN);
		
		for( auto spec : bioph.membraneProperties.capacitance_specs ) spec.apply(morph, segment_Cm_);
		for( auto spec : bioph.intracellularProperties.resistivity_specs ) spec.apply(morph, segment_Ra_);
		
		// TODO throw a mighty fit in case Cm and Ra are incomplete
		
		// 
		std::vector<Real> compartment_capacitance(number_of_compartments, NAN);
		
		// NB abuse NeuroML's tree morphology limitation to set axial resistance values in a compartment-parallel vector
		// where each value is child's axial resistance to its *parent*
		// Will need a more general way to store conductivity constants when
		//  non-tree neuron topologies are implemented LATER
		std::vector<Real> inter_compartment_axial_resistance(number_of_compartments, NAN);
		
		// Resistance is not derived by the total edge-to-edge resistance of compartments, 
		// because modellers ofter taper the ends of neurite to zero or another tiny value.
		// And this causes the half total resistance of a compartment to be disproportionally 
		// 	higher than the anticipated resistance of the proximal half of it,
		// 	since it accounts for the tapered end thich is not important for inter-compartment resistivity (since the neurite ends there).
		// Instead, the resistance between two compartments is the sum of the resistance 
		// 	of the distal half of the parent compartment and that of the proximal half of the child compartment.
		// Refer also to the code of NEURON that follows the same approach: treeset.cpp, routine diam_from_list
		// and to notes on the NEURON message board:
		// 	https://www.neuron.yale.edu/phpBB/viewtopic.php?f=15&t=2539&p=10078&hilit=axial+resistivity#p10078
		// 	https://www.neuron.yale.edu/phpBB/viewtopic.php?f=8&t=3904&p=16807&hilit=axial+resistivity#p16808
		// Note: This model is accurate for unbranched sections of neurite, but not for branches.
		// Conductance and area in that case are not that well defined, since integrating from the middle of the compartment to the precise fractionALong the compartment still may not include the small-scale details of the branch: a model should be made by an expert for that LATER, or more generally getting to override such conductance values. 
		// In any case, the true value of Ra 𓁛 is obscure anyway, so why should it matter... again, only an expert can tell.
		// Note: Compartments with zero length are assumed to have zero resistance, (as they have zero length) even though they have non-zero area and volume (as implied spheres).
		std::vector<Real> compartment_axial_resistance_proximal(number_of_compartments, NAN);
		std::vector<Real> compartment_axial_resistance_distal(number_of_compartments, NAN);
		
		for( int comp_seq = 0; comp_seq < number_of_compartments; comp_seq++ ){
			
			if(config.verbose){
				printf("Summating resistance and conductance for comp_seq %d\n", (int) comp_seq );
			}
			
			const auto comp_length = compartment_lengths[comp_seq];
			const auto comp_midpoint_length = comp_length / 2;
			
			auto &comp_axial_resistance_proximal = compartment_axial_resistance_proximal[comp_seq];
			auto &comp_axial_resistance_distal   = compartment_axial_resistance_distal  [comp_seq];
			compartment_axial_resistance_proximal[comp_seq] = 0;
			compartment_axial_resistance_distal  [comp_seq] = 0;

			auto &comp_capacitance = compartment_capacitance [comp_seq];
			comp_capacitance = 0;
			
			// also evaluate the d_lambda rule for the component, as a diagnostic for users who want to inspect the numerics at play
			const Real compartment_lambda_hertz = 100;
			Real compartment_lambda = 0;
			
			// summate by iterating through all segment-based fractions
			Real running_comp_length = 0;
			
			const auto &segs_in_comp = comp_disc.compartment_to_segment_seq[comp_seq];
			
			for( int seg_in_comp = 0; seg_in_comp < (int) segs_in_comp.size(); seg_in_comp++ ){
				Real from_fractionAlong = (seg_in_comp == 0) ? comp_disc.compartment_to_first_segment_start_fractionAlong[comp_seq] : 0;
				Real   to_fractionAlong = (seg_in_comp+1 == (int) segs_in_comp.size()) ? comp_disc.compartment_to_last_segment_end_fractionAlong[comp_seq] : 1;
				Int seg_seq = segs_in_comp[seg_in_comp];
				
				// TODO check if there is any case proximal should be the specified one, or the distal of parent seg instead
				
				if(config.verbose){
					printf("seg %d (id %d), fractionAlong from %g to %g\n", (int) seg_seq, (int) morph.lookupNmlId(seg_seq), from_fractionAlong, to_fractionAlong );
				}
				
				const Morphology::Segment &seg = morph.segments.atSeq(seg_seq);
				
				const Morphology::Segment::Point3DWithDiam seg_from = LerpPt3d(seg.proximal, seg.distal, from_fractionAlong);
				const Morphology::Segment::Point3DWithDiam seg_to   = LerpPt3d(seg.proximal, seg.distal,   to_fractionAlong);
				
				Real comp_fragment_length = DistancePt3d( seg_from, seg_to );
				Real comp_fragment_area   = GeomHelp::Area( comp_fragment_length, seg_from.d, seg_to.d );
				
				// Add capacitance to compartment fragment
				Real comp_fragment_capacitance = segment_Cm_[seg_seq] * comp_fragment_area;
				// and rescale to engine units
				comp_fragment_capacitance = ( (Scales<SpecificCapacitance>::native) * (microns*microns) ).ConvertTo( comp_fragment_capacitance, Scales<Capacitance>::native );
			
				comp_capacitance += comp_fragment_capacitance;
				
				// Add resistance to the proximal and distal halves of the compartment
				// running_comp_length is the sum of path length (arc length in NEURON lingo) of compartment fragments up to, excluding this iteration
				Real Ra = segment_Ra_[seg_seq];
				
				// input: Ra in native unit, start-end in microns, returns resistance in native unit
				auto GetFrustumResistance = [ &microns ]( Real Ra, Morphology::Segment::Point3DWithDiam from, Morphology::Segment::Point3DWithDiam to ){
					// using the volumetric frustum integral, instead of cylinder approximation
					// R = (Ra/pi)*( Length/(Radius_start*Radius_end) ) 
					// To get this formula, integrate dR = dx * Ra/CrossSection = dx * Ra /( pi * ( (1-x/Length)*Radius_start + (x/Length)*Radius_end )^2 ) for x = 0...Length
					Real resistance = Ra * (4/M_PI) * DistancePt3d(from, to) / ( from.d * to.d );
					// and rescale to engine units
					resistance = ( (Scales<Resistivity>::native * microns)/(microns^2) ).ConvertTo( resistance, Scales<Resistance>::native );
					return resistance;
				};
				
				if(config.verbose){
					printf("resistance measurement: from_comp_length %g to_comp_length %g midpoint_length %g\n", 
						running_comp_length, running_comp_length + comp_fragment_length, comp_midpoint_length
					);
				}
				// XXX this is a workaround for isolated zero length sections.
				// The advantages with this approach are that it:
				// - is simple to implement
				// - will behave like the NML exporter for the simple, expected cases (soma as sphere)
				// - can attempt to simulate any model, when the NML exporter's trick does not account for some contrived cases that will still get ill defined numerically. 
				// TODO follow the precise workaround that NeuroML does: If a segment has zero length AND proximal diameter = distal diameter AND ( it is the only segment in a cable  OR it is a self standing segment ) then consider it a cylinder with length = proximal diameter.
				// Accounting for path length on this cylinder should also have a small(insignificant?) effect on inhomogeneous distributions.
				// This also resolves the issue with whether Area should be based on a frustum or sphere.
				// TODO put warnings in place for improper tracing of neurons, to avoid halting when axial resistance predictably ends up invalid, when the specific conditions for the NML exporter's edge case do not apply (for example, if a segment has zero length but varying diameter and it is the only one in the cable, or when a cable includes multiple segments with zero diameter)
				// A variation could be to assume that all segments in a zero length cable are spheres...
				if( comp_length == 0 ){
					// Consider the form of a cylinder with diameter = length to calculate resistance.
					// Note that proximal and distal diameters amy differ, assume both are something sane...
					auto seg_to_extrapolated = seg_from;
					seg_to_extrapolated.y += seg_to.d;
					const Morphology::Segment::Point3DWithDiam seg_midpoint = LerpPt3d(seg.proximal, seg_to_extrapolated, 0.5);
					
					if(config.verbose){
						printf("\t zero length segment\n" );
					}
					
					comp_axial_resistance_proximal += GetFrustumResistance( Ra, seg_from, seg_midpoint);
					comp_axial_resistance_distal   += GetFrustumResistance( Ra, seg_midpoint, seg_to_extrapolated  );
				}
				else if( running_comp_length + comp_fragment_length < comp_midpoint_length ){
					// whole part of seg belongs to the proximal half
					if(config.verbose){
						printf("\t proximal half\n" );
					}
					comp_axial_resistance_proximal += GetFrustumResistance( Ra, seg_from, seg_to);
				}
				else if( comp_midpoint_length <= running_comp_length ){
					// whole part of seg belongs to the distal half
					if(config.verbose){
						printf("\t distal half\n" );
					}
					comp_axial_resistance_distal   += GetFrustumResistance( Ra, seg_from, seg_to);
				}
				else{
					assert( running_comp_length < comp_midpoint_length && comp_midpoint_length <= running_comp_length + comp_fragment_length );
					// seg straddles the proximal and distal half
					
					// NB: this is the fraction to lerp over the seg fraction, from seg_from to seg_to.
					// Not over the full extent of the seg, proximal to distal. 
					Real compartment_midpoint_fragment_fraction = (comp_midpoint_length - running_comp_length) / comp_fragment_length;
					// prefer NaN values to be mapped to fractionAlong = end, to keep the ilusion of zero length segments being processed whole
					if(!( compartment_midpoint_fragment_fraction < 1 )) compartment_midpoint_fragment_fraction = 1;
					if(!( 0 < compartment_midpoint_fragment_fraction )) compartment_midpoint_fragment_fraction = 0;
					
					const Morphology::Segment::Point3DWithDiam seg_midpoint = LerpPt3d(seg_from, seg_to, compartment_midpoint_fragment_fraction);
					
					if(config.verbose){
						printf("\t midpoint: fragment fraction %g\n", compartment_midpoint_fragment_fraction );
					}
					
					comp_axial_resistance_proximal += GetFrustumResistance( Ra, seg_from, seg_midpoint);
					comp_axial_resistance_distal   += GetFrustumResistance( Ra, seg_midpoint, seg_to  );
				}
				
				// Add the estimate of AC wave length as per the NEURON book, just as a diagnostic value. 
				
				// inputs in native units, except for time whoch is in Hz, and that is in order to omit it from the scaling factor... it seems to be wbe working nonetheless, why not leave it that way ... because Ra * Cm = Time/Length and that must be counterbalanced!!! The correctioun would be to set the ScaleEntry of Hz explicitly, and correct if it doesn't match Ra*Cm/microns XXX
				// Note that the average diameter is used in the calculation even though the integrand is a non linear function of diameter, if it's accurate enough for Hines it's good enough.
				auto lambda_f_microns = [&microns](float diam, float freq_Hz, float Ra, float Cm, bool debug){
					ScaleEntry scale_Ra = Scales<Resistivity>::native; // {"ohm_cm" ,-2, 1.0}; // NEURON units
					ScaleEntry scale_Cm = Scales<SpecificCapacitance>::native; // {"uF_per_cm2",-2, 1.0}; // NEURON units
					
					float lambda_f = std::sqrt( diam/( 4 * M_PI * freq_Hz * Ra * Cm ) );
					ScaleEntry dla_scale = ( microns / (scale_Ra * scale_Cm) )^(0.5);
					if( debug ){
						printf("dla %g %g %g\n", lambda_f, pow10(dla_scale.pow_of_10)*dla_scale.scale, lambda_f*(pow10(dla_scale.pow_of_10)*dla_scale.scale) );
					}
					return dla_scale.ConvertTo( lambda_f, microns );
				};
				
				Real comp_fragment_lambda_microns = lambda_f_microns( (seg_from.d + seg_to.d) / 2, compartment_lambda_hertz, Ra, segment_Cm_[seg_seq], config.debug );
				compartment_lambda += comp_fragment_length / comp_fragment_lambda_microns;
				
				// done for this compartment fragment, add path length to sum
				running_comp_length += comp_fragment_length;
			}
			
			if(config.verbose){
				printf("Compartment %d  capacitance: %g %s  proximal resistance: %g %s  distal resistance: %g %s\n", (int) comp_seq,
					comp_capacitance              , Scales<Capacitance>::native.name,
					comp_axial_resistance_proximal, Scales<Resistance >::native.name,
					comp_axial_resistance_distal  , Scales<Resistance >::native.name
				);
			}
			
			// report total 100Hz ac wavelength
			// try checking the d_lambda rule
			const float d_lambda = 0.1;
			// XXX why is this 0.9 factor here again ?? not a fixme since it matches the formula in NEURON
			float nseg_factor = ( compartment_lambda / d_lambda ) + 0.9;
			if( config.verbose ){
			printf("nseg %.9f\n", nseg_factor);
			}
			int nseg = int( nseg_factor / 2 ) * 2 + 1 ; //really should be odd, to avoid midpoint shenanigans
			(void) nseg; // It is useful as a diagnostic to show for modellers. subdivision perhaps LATER, but it would be better to use a separate neuron design tool for that
			
			// also take advantage of comp_seq order to resolve resistance with parent as well
			auto comp_parent = comp_disc.tree_parent_per_compartment[comp_seq];
			if(comp_parent >= 0){
				inter_compartment_axial_resistance[comp_seq] = comp_axial_resistance_proximal + compartment_axial_resistance_distal[comp_parent];
			
			}
			else{
				inter_compartment_axial_resistance[comp_seq] = NAN;
			}
		}
		if(config.verbose){
		for( int comp_seq = 0; comp_seq < number_of_compartments; comp_seq++ ){
			printf("Compartment %d resistance with parent %d: %g %s\n", 
				(int) comp_seq, (int) comp_disc.tree_parent_per_compartment[comp_seq],
				inter_compartment_axial_resistance[comp_seq] , Scales<Resistance>::native.name
			);
		}
		}
		for( int comp_seq = 0; comp_seq < number_of_compartments; comp_seq++ ){
			auto comp_parent = comp_disc.tree_parent_per_compartment[comp_seq];
			if( comp_parent < 0 ) continue;
			auto resistance = inter_compartment_axial_resistance[comp_seq];
			if(!( std::isfinite(resistance) && resistance > 0 )){
				// TODO prettier printing, though it shouldn't happen
				
				printf("internal error: Conductance between compartments %d, %d is undefined \n", 
				(int) comp_seq, (int) comp_parent );
				return false;
			};
		}
		
		// Approximate the smallest time constant of the passive system, using the Method of Time Constants. for RC circuits:
		// https://designers-guide.org/forum/Attachments/MTC.pdf
		// fastest pole = sum (pole of each capacitor with conductivities to ground, if all other capacitors were shorted)
		Real rate_total = 0;
		ScaleEntry RC_scale = (Scales<Resistance>::native * Scales<Capacitance>::native);
		for( int comp_seq = 0; comp_seq < number_of_compartments; comp_seq++ ){
			//conveniently, in neural compartment models, other capacitors being shorted means only directly adjacent resistances matter
			//  R1    R2    R3    R4
			// -vvv-+-vvv-+-vvv-+-vvv-
			//      |     =     |   
		    // GND -+-----+-----+- GND
			float Gtotal = 0;
			for( Int adjacent_comp : comp_connections[comp_seq] ){
				size_t R_index = (adjacent_comp > comp_seq) ? adjacent_comp : comp_seq ; // NB: remember how inter_compartment_axial_resistance is structured
				Real R = inter_compartment_axial_resistance[R_index];
				Gtotal += 1/R;
			}
			
			float rate = Gtotal / compartment_capacitance[comp_seq] ;
			rate_total += rate;
			
			float tau = RC_scale.ConvertTo( 1/rate, Scales<Time>::native); //TODO check values once more
			if( config.verbose ){
			printf(" compartment axial %g %s\n", tau, Scales<Time>::native.name );
			}
		}
		float tau_total = RC_scale.ConvertTo( 1/rate_total, Scales<Time>::native);
		if(config.verbose) printf(" total axial %g %s\n", tau_total, Scales<Time>::native.name );
		
		
		// Now perform analysis for Backward Euler method
		printf("\tAnalyzing Bwd Euler...\n");
		// Gauss elimination of one neuron's conductance matrix for Bwd Euler can be performed as a tree traversal, with the benefit of O(N) run time.
		// The process followed is thus equivalent to the Hines algorithm (1983).
		// The particular order of circuit nodes to visit could be tweaked. Here we use a depth-first traversal starting from node 0 (supposed to be the soma), but a BFS order could be used as well to reduce tree depth (and hopefully numerical error).
		struct BackwardEuler{
			public:
			// just a reference to the connection matrix to be traversed
			// typename could be made implicit with some template magic
			const std::vector< std::vector<int32_t> > &conn_list;
			// intermediate result: tha nodes visited by DFS.
			std::vector< bool > node_gray;
			// results: the sequence of nodes to visit, and the map of each node to its parent.
			std::vector< Int > order_list; // per processing step
			std::vector< Int > order_parent; // per node_seq
			// Recursivly traverse
			void DFS(
				Int i // node being visited
			){
				node_gray[i] = true;
				
				for( Int j : conn_list[i]){
					if( node_gray[j] ) continue;
					
					// otherwise, explore
					order_parent[j] = i;
					node_gray[j] = true;
					DFS( j );
				}
				
				order_list.push_back(i);
			}
			
			protected:
			BackwardEuler( const std::vector< std::vector<int32_t> > &_conn ) : conn_list(_conn) {
				
				Int nCells = (Int)conn_list.size();
				
				node_gray = std::vector< bool >( nCells, false );
				order_list = std::vector< Int >();
				order_parent = std::vector< Int >( nCells, -1 );
				
			}
			public:
			static void GetOrderLists(
				const std::vector< std::vector<int32_t> >conn_list,
				std::vector< Int > &order_list,
				std::vector< Int > &parent_list,
				Int start_from = 0
			){
				BackwardEuler ob(conn_list);
				
				ob.DFS( start_from );
				
				order_list = ob.order_list;
				parent_list = ob.order_parent;
			}
		};
		
		// Fill in the Bwd Euler order for the tree
		auto &order_list = pig.cable_solver.BwdEuler_OrderList;
		auto &order_parent = pig.cable_solver.BwdEuler_ParentList;
		BackwardEuler::GetOrderLists( comp_connections, order_list, order_parent );
		
		if( config.verbose ){
		printf("Order: ");
		for( Int val : order_list ){
			printf("%ld ", val);
		}
		printf("\n");
		printf("Parent: ");
		for( Int val : order_parent ){
			printf("%ld ", val);
		}
		printf("\n");
		}
		
		// Also fill in the part of the diagonal that's axial conductance to other, since it's going to be fixed during the simulation
		// 	(a made variable LATER)
		auto &compartment_adjacent_InvRC = pig.cable_solver.BwdEuler_InvRCDiagonal;
		compartment_adjacent_InvRC = std::vector<Real>(number_of_compartments, 0 );
		// and get the connectivity diagonals for bwd Euler
		for( int comp_seq = 0; comp_seq < number_of_compartments; comp_seq++ ){
			for( int adjacent_comp : comp_connections[comp_seq] ){
				Int idx = std::max( comp_seq, adjacent_comp ); // NB: remember how inter_compartment_axial_resistance is structured
				auto R = inter_compartment_axial_resistance[ idx ];
				auto C = compartment_capacitance[ comp_seq ];
				auto D = ( (Scales<Resistance>::native * Scales<Capacitance>::native)^(-1) ).ConvertTo( 1/(R*C), Scales<Frequency>::native );
				compartment_adjacent_InvRC[comp_seq] += D;
			}
		}
		
		if( config.verbose ){
		printf("Diagonal 1/RC Constant(%s): ", Scales<Frequency>::native.name );
		for( auto val : compartment_adjacent_InvRC ){
			printf("%g ", val);
		}
		printf("\n");
		}
		
		// initial potential, threshold for every compartment
		std::vector<float> compartment_V0(number_of_compartments, NAN);
		std::vector<float> compartment_Vt(number_of_compartments, NAN); // TODO verify it is not NAN when a spiking connection is defined in respective compartment
		// TODO throw a mighty fit in case V0 is incomplete
		
		// the mapping from segments to compartments is not necessarily direct.
		// 	Produce a list that contains all compartments that overlap with said segments
		auto SpecToCompList = [](
			const Morphology &morph, const AcrossSegOrSegGroup &spec, 
			const CellInternalSignature::PhysicalCell::CompartmentDiscretization &comp_disc
		){
			IdListRle comp_list;
			spec.reduce( morph, [ &comp_disc, &comp_list ]( Int seg_seq ){
				for( Int comp_seq : comp_disc.segment_to_compartment_seq[seg_seq] ) comp_list.Addd(comp_seq);
			} );
			comp_list.Compact();
			return comp_list;
		};
		auto ApplyFixedValue = []( const IdListRle &id_list, const auto &value, auto &array ){
			id_list.reduce( [&value, &array]( Int comp_seq ){
				array[comp_seq] = value;
			} );
		};
		
		for( auto spec : bioph.membraneProperties.initvolt_specs ) ApplyFixedValue(SpecToCompList(morph, spec, comp_disc), spec.value, compartment_V0);
		for( auto spec : bioph.membraneProperties.threshold_specs ) ApplyFixedValue(SpecToCompList(morph, spec, comp_disc), spec.value, compartment_Vt);
		
		// Realize per-compartment signatures
		
		auto &comp_definitions = pig.comp_definitions;
		comp_definitions.resize(number_of_compartments);
		
		for( int comp_seq = 0; comp_seq < number_of_compartments; comp_seq++ ){
			auto &comp_def = comp_definitions[comp_seq]; 
			comp_def.V0 = compartment_V0[comp_seq];
			comp_def.Vt = compartment_Vt[comp_seq];
			comp_def.AxialResistance = inter_compartment_axial_resistance[comp_seq];
			comp_def.Capacitance = compartment_capacitance[comp_seq];
			
			comp_def.adjacent_compartments = comp_connections[comp_seq];
		}
		
		// add the ion channel distributions too
		for(const auto &spec : bioph.membraneProperties.channel_specs){
			
			auto comp_arr = SpecToCompList(morph, spec, comp_disc).toArray();
			
			// used when the is not uniform uniform per compartment
			std::vector<Real> inho_parameter_values;
			
			if(spec.conductivity.type == ChannelDistribution::Conductivity::NON_UNIFORM){
				inho_parameter_values.assign(comp_arr.size(), NAN);
				
				assert( ((const AcrossSegOrSegGroup &) spec).type == AcrossSegOrSegGroup::Type::GROUP );
				const auto &seg_group = morph.segment_groups[spec.seqid];
				
				const auto &inho = spec.conductivity.inho;
				assert( inho.parm >= 0 );
				const auto &inhoparm = seg_group.inhomogeneous_parameters.get(inho.parm);
				
				// Fetch the parameter values
				if( inhoparm.metric == Morphology::InhomogeneousParameter::PATH_LENGTH_FROM_ROOT ){
					for( int32_t comp_in_spec = 0; comp_in_spec < (int32_t)comp_arr.size(); comp_in_spec++ ){
						int32_t comp_seq = comp_arr[comp_in_spec];
						inho_parameter_values[comp_in_spec] = compartment_path_length_from_root[comp_seq];
					}
				}
				else{
					printf("internal error: unknown inhomogeneous parameter type %d", (int)inhoparm.metric );
					return false;
				}
				
				// Optionally post-process in the context of this segment group's set of values.
				// p0 is the minimum value found in the set of samples, and p1 is the maximum found. 
				// These "guard" values are chosen intentionally, to match NEURON's SubsetDomainIterator behaviour. (Refer to subiter.hoc file)
				Real p0 = 1e9;
				Real p1 = 1;
				for( auto &val : inho_parameter_values ) p0 = std::min(p0, val);
				for( auto &val : inho_parameter_values ) p1 = std::max(p1, val);
				
				if( inhoparm.subtract_the_minimum ){
					for( auto &val : inho_parameter_values ) val -= p0;
					p1 -= p0;
					p0 = 0;
				}
				if( inhoparm.divide_by_maximum ){
					for( auto &val : inho_parameter_values ) val /= p1;
					p0 /= p1;
					p1 = 1;
				}
				
				// Also make sure that random() is not used in the arithmetic expression of the distribution
				for( const auto &term : inho.value.terms ){
					if( term.type == Term::Type::RANDOM ){
						printf("inhomogeneous ion channel conductivity involving random() not supported yet\n");
						return false;
					}
				}
			}
			
			for( int32_t comp_in_spec = 0; comp_in_spec < (int32_t)comp_arr.size(); comp_in_spec++ ){
				int32_t comp_seq = comp_arr[comp_in_spec];
				
				CellInternalSignature::IonChannelDistributionInstance instance;
				
				instance.ion_species = spec.ion_species;
				instance.ion_channel = spec.ion_channel;
				instance.type = spec.type;
				
				if(spec.conductivity.type == ChannelDistribution::Conductivity::FIXED){
					instance.conductivity = spec.conductivity.value;
				}
				else if(spec.conductivity.type == ChannelDistribution::Conductivity::NON_UNIFORM){
					// XXX the inhomogeneous parameters may contain random(). Thus, in the general case, they have to be evaluated separately for each cell, at instantiation time !
					
					const auto &inho = spec.conductivity.inho;
					
					// the only symbol is the inhomogeneous parameter, and it could be either dimensionless microns or dimensionless pure, 
					// depending on if inhoparm.divide_by_maximum was used 
					Real value_context[1] = { inho_parameter_values[comp_in_spec] };
					const Dimension dimension_context[1] = { Dimension::Unity() };
					assert( inho.value.symbol_refs.size() <= 1 );
					
					Real val = NAN;
					Dimension dim_unused = Dimension::Unity(); // should be dimensionless SI conductivity when parsed
					Int random_call_counter = 0; // no other calls to random() that need to be put on a counter
					EvaluateLemsExpression(
						inho.value, value_context, dimension_context,
						dimensions, random_call_counter, val, dim_unused, 
						[]( Real random_arg, Int random_call_counter ){
							(void) random_arg; (void) random_call_counter;
							// random numbers might be used LATER, only if necessary. In this scope, we can't really use a value because Evaluate must be called seperately for each neuron instance, at instantiation stage ... not codegen stage
							return NAN;
						}
					);
					
					// val is given in SI value, convert it to engine value
					const ScaleEntry mhos_per_m2 = {"S_per_m2"  , 0, 1.0};
					// assume that the returned value is always SI conductivity whether the factor is length or normalised, TODO explain this somewhere!
					// printf("evall %10g %10g %10g\n",value_context[0],val, mhos_per_m2.ConvertTo( val, Scales<Conductivity>::native) );
					
					instance.conductivity = mhos_per_m2.ConvertTo( val, Scales<Conductivity>::native );
				}
				else{
					printf("internal error: unknown inhomogeneous ion channel conductivity type\n");
					return false;
				}
				
				instance.erev = spec.erev;
				instance.vshift = spec.vshift;
				instance.permeability = spec.permeability;
				instance.number = spec.number;
				
				
				comp_definitions[comp_seq].ionchans.push_back(instance);
			};
		}
		
		// and the concentration models
		for(const auto &spec : bioph.intracellularProperties.ion_species_specs){
			CellInternalSignature::IonSpeciesDistributionInstance instance;
			
			//instance.ion_species = spec.species;
			instance.conc_model_seq = spec.concentrationModel;
			
			instance.initialConcentration = spec.initialConcentration;
			instance.initialExtConcentration = spec.initialExtConcentration;
			
			SpecToCompList(morph, spec, comp_disc).reduce([&comp_definitions, &spec, instance]( Int comp_seq ){
				comp_definitions[comp_seq].ions[spec.species] = instance;
			});
		}
		// what TODO with external concentrations ??
		
		// fill in input/synapse occurences to compartment signatures
		for( const auto &keyval : input_types_per_cell_per_compartment[cell_seq] ){
			comp_definitions[keyval.first].input_types = keyval.second;
		}
		for( const auto &keyval : synaptic_component_types_per_cell_per_compartment[cell_seq] ){
			comp_definitions[keyval.first].synaptic_component_types = keyval.second;
		}
		for( const auto &key : spiking_outputs_per_cell_per_compartment[cell_seq] ){
			comp_definitions[key].spike_output = true;
		}
		
		
		// pick a cable equation integrator
		
		auto &cell_cable_solver = pig.cable_solver.type = config.cable_solver;
		if( cell_cable_solver == SimulatorConfig::CABLE_SOLVER_AUTO ){
			cell_cable_solver = SimulatorConfig::CABLE_BWD_EULER; // TODO
		}
		// bool postupdate_inside_cell = false; // LATER
		
		// cell analysis complete
		
		// now compose variables for the work unit
		
		// split the signature-construction part here, to separate per-cell from per-compartment analysis

		
		// constants
		size_t Index_Capacitance = AppendSingle_CellScope.Constant( compartment_capacitance, "Compartment Capacitance ("+std::string(Scales<Capacitance>::native.name)+")" );
		size_t Index_AxialResistance = AppendSingle_CellScope.Constant( inter_compartment_axial_resistance, "Axial Resistance ("+std::string(Scales<Resistance>::native.name)+")" );
		size_t Index_VoltageThreshold = AppendSingle_CellScope.Constant( compartment_Vt, "Spike Threshold ("+std::string(Scales<Voltage>::native.name)+")" );
		size_t Index_MembraneArea = AppendSingle_CellScope.Constant( compartment_areas, "Membrane Surface Area (microns^2)" );
		size_t Index_Temperature = AppendSingle_CellScope.Constant( net.temperature, "Temperature (K)" ); // TODO move to global
		
		// and states
		pig.Index_Voltages = AppendSingle_CellScope.StateVariable( compartment_V0, "Voltage ("+std::string(Scales<Voltage>::native.name)+")"  );
		
		// more constants and variables will arise, as the code is composed
		
		// const std::string Rate_suffix = Convert::Suffix(
		// 	( Scales<Frequency>::native * Scales<Time>::native ) // to unitless
		// );
		
		// now on to parts of the physical cell
		
		printf("Generating code for %s...:\n", sig.name.c_str());
		char tmps[10000]; // buffer for a single code line
		
		EmitKernelFileHeader( sig.code );
		EmitWorkItemRoutineHeader( sig.code );
		
		
		const std::string tab = "\t";
		
		// some per-cell stuff
		
		sig.code += CloneSubitemIndices("cell", "local", "\t");
		sig.code += ExposeSubitemContext("cell", "global", "\t");
		
		sig.code += "	\n";
		sprintf(tmps, "	const float temperature = cell_constants[%zd]; //a global if there ever was one\n", Index_Temperature); sig.code += tmps;
		
		sig.code +=   "	\n";
		sprintf(tmps, "	const float *V = &cell_state[%zd]; \n", pig.Index_Voltages); sig.code += tmps;
		sprintf(tmps, "	      float *V_next = &cell_stateNext[%zd]; \n", pig.Index_Voltages); sig.code += tmps;
		sprintf(tmps, "	const float *R_Axial = &cell_constants[%zd]; \n", Index_AxialResistance); sig.code += tmps;
		sprintf(tmps, "	const float *C = &cell_constants[%zd]; \n", Index_Capacitance); sig.code += tmps;
		sprintf(tmps, "	const float *V_threshold = &cell_constants[%zd]; \n", Index_VoltageThreshold); sig.code += tmps;
		sprintf(tmps, "	const float *Area = &cell_constants[%zd]; \n", Index_MembraneArea); sig.code += tmps;
		
		// also append diagonal conductance matrix entries for the Hines algorithm
		if( cell_cable_solver == SimulatorConfig::CABLE_BWD_EULER ){
		auto &cabl_impl = pig.cable_solver_implementation;
		cabl_impl.Index_BwdEuler_InvRCDiagonal = AppendMulti_CellScope.Constant("Bwd Euler Diagonal 1/RC Constant");
		cabl_impl.Index_BwdEuler_WorkDiagonal = AppendMulti_CellScope.StateVariable("Bwd Euler Diagonal Scratchpad");
		sprintf(tmps, "	const Table_F32 PerComp_InvRC_Axial  = cell_const_table_f32_arrays[%zd];\n", cabl_impl.Index_BwdEuler_InvRCDiagonal); sig.code += tmps;
		sprintf(tmps, "	Table_F32 PerComp_InvRC_Diagonal = cell_state_table_f32_arrays[%zd];\n", cabl_impl.Index_BwdEuler_WorkDiagonal); sig.code += tmps;
		}
		
		sig.code += "	\n";
		
		
		// per cell RNG
		ImplementRngSeed(
			AppendSingle_CellScope,
			"", tab,
			"cell",
			sig.common_in_cell.cell_rng_seed,
			sig.code
		);
		sig.code += "	const int rng_object_id = cell_rng_seed;\n";
		
		sig.code += "	\n";
		
		//open up this table too
		pig.comp_implementations.resize( number_of_compartments );
		
		// first, evaluate the internal chemical dynamics of each compartment, and integrate them inline
		auto ImplementInternalCompartmentIntegration = [
				&config,
				&model, &ion_channels, &conc_models, &ion_species, &dimensions, &component_types, &microns,
				&ImplementSynapseType, &ImplementInputSource
		](
			const SignatureAppender_Single &AppendSingle, const SignatureAppender_Table &AppendMulti,
			const InlineLems_AllocatorCoder &DescribeLemsInline,
			const std::string &for_what,
			const std::string &tab,
			bool  flatten_adjacency,
			const SimulatorConfig::CableEquationSolver &cell_cable_solver,
			const BiophysicalProperties &bioph,
			const CellInternalSignature::PhysicalCell::CompartmentDefinition &comp_def,
			CellInternalSignature::PhysicalCell::CompartmentImplementation &comp_impl,
			Int &random_call_counter,
			std::string &ccde
		){
			char tmps[10000];
			auto AppendConstant = [&AppendSingle]( Real default_value, const std::string &for_what){
				return AppendSingle.Constant( default_value, for_what );
			};
			auto AppendStateVariable = [&AppendSingle]( Real default_value, const std::string &for_what){
				return AppendSingle.StateVariable( default_value, for_what );
			};
			
			// ccde += "\t{\n";
			
			sprintf(tmps, "	float Acomp = Area[comp];\n"); ccde += tmps; // in microns^3 TODO verify once more
			
			
			sprintf(tmps, "	float Vcomp = V[comp];\n"); ccde += tmps;
			
			bool uses_Iaxial = (
				cell_cable_solver == SimulatorConfig::CABLE_FWD_EULER
			);
			
			// if such quantities must be communicated, they should be put in a scratch vector for compartment currents
			sprintf(tmps, "	float I_internal = 0;\n"); ccde += tmps; // total inward current of membrane mehcanisms (ie not axial)
			// for now, integrate voltage on the spot (fwd euler can be converted to bwd euler this way just by correcting by + Gtotal * V) to store this per-compartment value
			
			sprintf(tmps, "	float G_internal = 0;\n"); ccde += tmps; // total conductance of internal mechanisms, coupling to GND
			// Sum of the Thévenin equivalent series conductance of each one-port mechanism connecting the membrane to ground (such as ion channels, post-synapses, gap junctions, conductance clamps...)
			// NOTE though that sorts of coupling that are not handled by a specific solver may also be binned to this, as Vext and associated are is assumed fixed (e.g. for gap junctions possibly)
			// for now, use the scratchpad diagonal vector to store this per-compartment value
			
			// TODO also use a 'total conductivity' variable to avoid the inefficiency and precision loss of multiplying by Area before summing these factors
			
			
			// peek into ion species to get some species populations
			// it's actually more elegant to contribute to ions, than to poll ion sources
			(void) ion_species; // LATER this might be useful
			ccde += "	// Ion flux sources\n";
			for( auto keyval : comp_def.ions){
				sprintf(tmps, "		float I_ion_%ld = 0; //total ion current\n", keyval.first); ccde += tmps;
				sprintf(tmps, "		float Conc_ion_%ld_intra = 0; //ion concentration intra\n", keyval.first); ccde += tmps;
				sprintf(tmps, "		float Conc_ion_%ld_extra = 0; //ion concentration extra\n", keyval.first); ccde += tmps;
			}
			
			auto ExposeRequirements_ConcModel = [  ]( Int ion_seq, const auto &distimpl, const std::string &tab){
				// internally expose stuff that is required for LEMS components
				std::string ret;
				ret += tab + "float Iion = I_ion_"+itos(ion_seq)+";\n";
				ret += tab + "float InitConcIntra = local_constants["+itos(distimpl.Index_InitIntra)+"];\n";
				ret += tab + "float InitConcExtra = local_constants["+itos(distimpl.Index_InitExtra)+"];\n";
				return ret;
			};
			
			// Ion concentrations should also be defined here
			// NB concentration exposure should not be a function of ion flux, or a circular dependency will form between concentration and ion channel components !
			ccde += "	// Ion concentrations\n";
			for( auto keyval : comp_def.ions){
				
				Int species_seq = keyval.first;
				const auto &instance = keyval.second;
				
				// now work with ion pools
				CellInternalSignature::IonSpeciesDistImplementation distimpl;
				
				const auto &conc_model = conc_models.get(instance.conc_model_seq);
				
				const std::string &tab = "\t";
				std::string ionpool_code;
				char tmps[2000];
				const std::string &for_that = for_what;
				std::string for_what = for_that+" Ion "+itos(species_seq)+" pool";
				
				
				// allocate some constants for every distribution ie. initial states
				
				distimpl.Index_InitIntra = AppendConstant(instance.initialConcentration, for_what+" Initial Internal Concentration");
				distimpl.Index_InitExtra = AppendConstant(instance.initialExtConcentration, for_what+" Initial External Concentration");
				
				ionpool_code += tab + "{\n";
				ionpool_code += tab;
				
				// TODO deduplicate also in dynamics section
				ionpool_code += ExposeRequirements_ConcModel( species_seq, distimpl, tab );
				if(conc_model.type == ConcentrationModel::COMPONENT){
					ionpool_code += "// LEMS component\n";
					const auto &comptype = model.component_types.get(conc_model.component.id_seq);
					
					distimpl.component = DescribeLems::AllocateSignature(comptype, conc_model.component, &AppendSingle, for_what + " LEMS");
					std::string lemscode = DescribeLems::Assigned(comptype, dimensions, distimpl.component, &AppendSingle, for_what, tab, random_call_counter);
					ionpool_code += lemscode;
					
					// update later on
					// expose whatever's possible here
					ionpool_code += DescribeLems::Exposures(comptype, for_what, tab, config.debug);
					
					// expose the exposures
					sprintf(tmps, "Conc_ion_%ld_intra = Lems_exposure_concentration;\n", species_seq); ionpool_code += tab+tmps;
					sprintf(tmps, "Conc_ion_%ld_extra = Lems_exposure_extConcentration;\n", species_seq); ionpool_code += tab+tmps;
				}
				else{
					// it is a built-in type
					
					// constants
					distimpl.Index_RestConc = AppendConstant(conc_model.restingConc, for_what+" Resting Concentration" );
					distimpl.Index_DecayTau = AppendConstant(conc_model.decayConstant, for_what+" Decay Tau" );
					std::string leak_factor_name = "LeakFactor???";
					if(conc_model.type == ConcentrationModel::LEAKY){
						leak_factor_name = "Shell Thickness";
					}else if(conc_model.type == ConcentrationModel::FIXED_FACTOR) {
						leak_factor_name = "Rho Factor";
					}
					distimpl.Index_Shellthickness_Or_RhoFactor = AppendConstant(conc_model.shellThickness_or_rhoFactor, for_what + " " + leak_factor_name );
					
					
					// and states
					distimpl.Index_Intra = AppendStateVariable(instance.initialConcentration   , for_what + "Intra" );
					// printf("index intra %zd\n", distimpl.Index_Intra);
					distimpl.Index_Extra = AppendStateVariable(instance.initialExtConcentration, for_what + "Extra" );
					
					// exposures are straightforward
					sprintf(tmps, "Conc_ion_%ld_intra = local_state[%zd];\n", species_seq, distimpl.Index_Intra); ionpool_code += tab+tmps;
					sprintf(tmps, "Conc_ion_%ld_extra = local_state[%zd];\n", species_seq, distimpl.Index_Extra); ionpool_code += tab+tmps;
					
				}
				ionpool_code += tab + "}\n";
				ccde += ionpool_code;
				
				comp_impl.concentration[species_seq] = distimpl;
				
			}
			
			// also check for the blessed calcium concentrations and fluxes
			auto Ca_species_seq  = bioph.Ca_species_seq ;
			auto Ca2_species_seq = bioph.Ca2_species_seq;
			//printf("\n\n\n\ncalcium %ld\n", Ca_species_seq);
			
			// LATER add a variable set of ions in a flexible way, I guess?
			bool has_Ca = false, has_Ca2 = false;
			if( comp_def.ions.count(Ca_species_seq) ){
				sprintf(tmps, "	const float Ca_concentration = Conc_ion_%ld_intra;\n", Ca_species_seq ); ccde += tmps;
				sprintf(tmps, "	const float Ca_concentration_extra = Conc_ion_%ld_extra;\n", Ca_species_seq ); ccde += tmps;
				has_Ca = true;
			}
			else{
				// omit
				sprintf(tmps, "	const float Ca_concentration = 0;\n"); ccde += tmps;
				sprintf(tmps, "	const float Ca_concentration_extra = 0;\n"); ccde += tmps;
				has_Ca = false;
			}
			if( comp_def.ions.count(Ca2_species_seq) ){
				sprintf(tmps, "	const float Ca2_concentration = Conc_ion_%ld_intra;\n", Ca2_species_seq ); ccde += tmps;
				sprintf(tmps, "	const float Ca2_concentration_extra = Conc_ion_%ld_extra;\n", Ca2_species_seq ); ccde += tmps;
				has_Ca = true;
			}
			else{
				// omit
				sprintf(tmps, "	const float Ca2_concentration = 0;\n"); ccde += tmps;
				sprintf(tmps, "	const float Ca2_concentration_extra = 0;\n"); ccde += tmps;
				has_Ca = false;
			}
			
			
			if( !has_Ca ){
				// could complain, but it's best to assume zero concentration for compatibility with jNeuroML
			}
			if( !has_Ca2 ){
				// could complain, but it's best to assume zero concentration for compatibility with jNeuroML
			}
			
			// TODO the cable equation handling should be broken off since it involves other compartments, therefore it is not really internal (important for compartment as work itam)
			if( uses_Iaxial ){
				const std::string Iaxial_suffix = Convert::Suffix(
					(Scales<Voltage>::native / Scales<Resistance>::native)
					.to(Scales<Current>::native)
				);
				
				ccde += tab+"// Inter-compartment leaks\n";
				ccde += tab+"	float I_axial = 0;\n";
				ccde += tab+"	int adj_conductance = -1;\n";
				ccde += tab+"	int adj_comp = -1;\n";
				
				std::string Adjcon_line = tab+"if( adj_conductance < comp ) adj_conductance = comp;\n";
				std::string Iaxial_line = tab+"I_axial += ( (V[adj_comp] - V[comp]) / R_Axial[adj_conductance] )"+Iaxial_suffix+";\n";
				std::string Both_lines = Adjcon_line + Iaxial_line;
				
				if( flatten_adjacency ){
					// add leaks (and perhaps longitudinal diffusions LATER)
					
					for( Int adjacent_comp : comp_def.adjacent_compartments ){
						ccde += tab+"adj_comp = "+itos(adjacent_comp)+"; \n";
						ccde += tab+"adj_conductance = adj_comp; \n";
						
						ccde += Both_lines;
					}
					ccde += tab+"// adj_conductance conditional should be optimized out in flattened code\n";
				}
				else{
					// also necessary for compartment deduplication
					auto Index_AdjComp = comp_impl.Index_AdjComp = AppendMulti.ConstI64( for_what+"Adjacent Compartments" );
					ccde += tab+ "const Table_I64 AdjCompartments = local_const_table_i64_arrays["+itos(Index_AdjComp)+"];\n";
					ccde += tab+ "const long long AdjComp_Count = local_const_table_i64_sizes["+itos(Index_AdjComp)+"];\n";
					ccde += tab+ "for( long long adjcomp_idx = 0; adjcomp_idx < AdjComp_Count; adjcomp_idx++ ){\n";
					
					ccde += tab+"\t""int adj_comp = AdjCompartments[adjcomp_idx];\n";
					ccde += tab+"\t""int adj_conductance = adj_comp; \n";
					
					ccde += Both_lines;
					
					ccde += tab+ "}\n";
				}
			}
			
			// TODO add ion channel, synapse etc. conductances to backward solver's conductivity
			
			//add ion channels
			ccde += "	// Current from ion channels\n";
			sprintf(tmps, "	float I_channels_total = 0;\n"); ccde += tmps;
			sprintf(tmps, "	float G_channels_total = 0;\n"); ccde += tmps;
			
			comp_impl.channel.resize(comp_def.ionchans.size());
			
			for( size_t inst_seq = 0; inst_seq < comp_def.ionchans.size(); inst_seq++ ){
				const auto &inst = comp_def.ionchans[inst_seq];
				const IonChannel &chan = ion_channels.get(inst.ion_channel);
				
				auto &distimpl = comp_impl.channel[inst_seq];
				
				const auto &for_that = for_what;
				std::string for_what = for_that +" ChannelDist "+itos(inst_seq);
					
				ccde   += "	{\n";
				
				// the only thing distribution types differ in is Erev, and possibly passing a VShift parameter to gates 
				// vShift should be available and set to zero, for compatibility, as the original LEMS types prescribe
				// TODO stop simulations from running with an unset (zero) caclium concentration, in that case
				
				
				// Fixed    uses a fixed Gbase and uses a fixed Erev
				// vShift is just like Fixed
				// Nernst   uses a fixed Gbase and provides own Erev
				// GHK1          uses no Gbase and provides own current
				// GHK2     uses a fixed Gbase and provides own current
				// Population   gets ohm Gbase and uses a fixed Erev
				bool is_population = ( inst.type == ChannelDistribution::POPULATION );
				
				bool uses_conductivity = (
					inst.type == ChannelDistribution::FIXED
					|| inst.type == ChannelDistribution::VSHIFT
					|| inst.type == ChannelDistribution::NERNST
					|| inst.type == ChannelDistribution::NERNST_CA2
					|| inst.type == ChannelDistribution::GHK2
				);
				
				// uses Gbase somewhere, even if modifying it
				bool uses_fixed_conductivity = uses_conductivity;
				
				bool provides_current = is_population;
				
				bool provides_density = !provides_current;
				
				bool fixed_Erev = (
					inst.type == ChannelDistribution::FIXED
					|| inst.type == ChannelDistribution::VSHIFT
					|| inst.type == ChannelDistribution::POPULATION
				);
				
				// a first pass through channel distributions, to initialize values such as vShift
				// ccde   += "	float Erev = NAN;\n";
				ccde   += "	float Vshift = 0;\n";
				
				if( inst.type == ChannelDistribution::VSHIFT ){
					std::size_t Index_Channel_Vshift = AppendConstant( inst.vshift, for_what + " Vshift" );
					sprintf(tmps, "		Vshift  = local_constants[%zd];\n", Index_Channel_Vshift); ccde += tmps;
					// FIXME implement Vshift
					printf("internal error: Vshift not yet implemented\n");
				}
				
				if( fixed_Erev ){
					//add constants, well, to Constants
					std::size_t Index_Channel_Erev = AppendConstant( inst.erev, "Erev for Fixed channel "+itos(inst_seq) );
					sprintf(tmps, "		float Erev  = local_constants[%zd];\n", Index_Channel_Erev); ccde += tmps;
					
				}
				else if(
					inst.type == ChannelDistribution::NERNST
					|| inst.type == ChannelDistribution::NERNST_CA2
				){
					double R = 8.3144621; // ( J / ( K * mol ) )
					double zCa = 2; // any ion LATER
					double F = 96485.3; // ( Cb / mol )
					auto SI_To_Erev_suffix = Convert::Suffix( Scales<Dimensionless>::native.to( Scales<Voltage>::native ) ); // SI (fundamental, really) units to Erev
					if( inst.type == ChannelDistribution::NERNST_CA2 ){
						sprintf(tmps, "		float Erev  = ( (float)(%.17g) * temperature / (float)( %.17g * %.17g) * logf( Ca2_concentration_extra / Ca2_concentration )%s );\n", R, zCa, F, SI_To_Erev_suffix.c_str()); ccde += tmps;
					}
					else{
						sprintf(tmps, "		float Erev  = ( (float)(%.17g) * temperature / (float)( %.17g * %.17g) * logf( Ca_concentration_extra / Ca_concentration )%s );\n", R, zCa, F, SI_To_Erev_suffix.c_str()); ccde += tmps;
					}
				}
				else if(
					inst.type == ChannelDistribution::GHK
					|| inst.type == ChannelDistribution::GHK2
				){
					// no Erev
				}
				else{
					printf("internal error: ion channel distribution not specifying use of Erev %d\n", inst.type);
					return false;
				} // end switch native ion channel type
				
				
				ccde   += "	float ChannelOpenFraction = NAN;\n";
				ccde   += "	float ChannelConductance = NAN;\n";
				// implement ion channel
				if( chan.type == IonChannel::COMPONENT ){
					ccde   += "	{\n";
					ccde   += DescribeLemsInline.SingleInstance( chan.component, "\t", for_what, distimpl.channel_component, config.debug );
					sprintf(tmps, "	ChannelOpenFraction = Lems_exposure_fcond;\n"); ccde += tmps;
					if( component_types.get(chan.component.id_seq).common_exposures.conductance >= 0 ){
						sprintf(tmps, "	ChannelConductance = Lems_exposure_g;\n"); ccde += tmps;
					}
					// TODO add support for conductivity
					ccde   += "	}\n";
				}
				else{
					// build a native channel
					distimpl.per_gate.resize(chan.gates.contents.size());
					
					struct DescribeRateThing{
						
						static double Value(float voltage_engine_units, const IonChannel::Rate &rate){
							const auto &Vcomp = voltage_engine_units;
							if(rate.type == IonChannel::Rate::EXPONENTIAL){
								return rate.formula.rate * exp( (Vcomp - rate.formula.midpoint ) / rate.formula.scale ) ; //TODO suffix
							}
							else if(rate.type == IonChannel::Rate::EXPLINEAR){
								double x = (Vcomp - rate.formula.midpoint) / rate.formula.scale;
								if(x == 0) return rate.formula.rate;
								else return rate.formula.rate * x / ( 1 - exp( -x ) );
							}
							else if(rate.type == IonChannel::Rate::SIGMOID){
								return  rate.formula.rate / (1 + exp( (rate.formula.midpoint - Vcomp ) / rate.formula.scale ) );
							}
							else{
								return NAN; // TODO handle with LEMS or sth?
							}
						}
					};
					auto DescribeRate_Thing = [&]( const IonChannel::Rate &rate, const std::string &tab, const std::string &for_what, const char *thing_name,  CellInternalSignature::ComponentSubSignature &component ){
						std::string rate_code;
						char tmps[2000];
						
						rate_code += tab; rate_code += "float "+std::string(thing_name)+"; // define exposure\n";
							
						if(rate.type == IonChannel::Rate::COMPONENT){
							rate_code += DescribeLemsInline.SingleInstance( rate.component, tab, for_what, component, config.debug );
							sprintf(tmps, "%s = Lems_exposure_%s;\n", thing_name, thing_name); rate_code += tab+tmps;
						}
						else{
							// it is a built-in type
							sprintf(tmps, "%s = ", thing_name); rate_code += tab+tmps;
							if(
								rate.type == IonChannel::Rate::EXPONENTIAL
								|| rate.type == IonChannel::Rate::EXPLINEAR
								|| rate.type == IonChannel::Rate::SIGMOID
							){
								std::size_t Index_Gate_BaseRate = AppendConstant(rate.formula.rate,     for_what + " Base" );
								std::size_t Index_Gate_Midpoint = AppendConstant(rate.formula.midpoint, for_what + " Mid" );
								std::size_t Index_Gate_Scale    = AppendConstant(rate.formula.scale,    for_what + " Scale");
							
								if(rate.type == IonChannel::Rate::EXPONENTIAL){
									sprintf(tmps, "local_constants[%zd] * expf( (Vcomp - local_constants[%zd] ) / local_constants[%zd] );\n", Index_Gate_BaseRate, Index_Gate_Midpoint, Index_Gate_Scale); rate_code += tmps;
								}
								else if(rate.type == IonChannel::Rate::EXPLINEAR){
									sprintf(tmps, "local_constants[%zd] * ( ( Vcomp == local_constants[%zd]) ? 1 : ( ( (Vcomp - local_constants[%zd] ) / local_constants[%zd] )  / (1 - expf( - (Vcomp - local_constants[%zd] ) / local_constants[%zd] ) ) ) );\n", Index_Gate_BaseRate, Index_Gate_Midpoint, Index_Gate_Midpoint, Index_Gate_Scale, Index_Gate_Midpoint, Index_Gate_Scale); rate_code += tmps;
								}
								else if(rate.type == IonChannel::Rate::SIGMOID){
									sprintf(tmps, "local_constants[%zd] / (1 + expf( (local_constants[%zd] - Vcomp ) / local_constants[%zd] ) );\n", Index_Gate_BaseRate, Index_Gate_Midpoint, Index_Gate_Scale); rate_code += tmps;
								}
							}
							else if( rate.type == IonChannel::Rate::FIXED ){
								std::size_t Index_Gate_Constant = AppendConstant(rate.formula.constant, for_what + " Fixed");
								sprintf(tmps, "local_constants[%zd];\n", Index_Gate_Constant); rate_code += tmps;
							}
							else{
								printf("internal error: ion channel rate thing type %s\n", itos(rate.type).c_str());
								assert(false);
							}
						}
						return rate_code;
					};
					
					auto DescribeRate_Rate = [&DescribeRate_Thing](const IonChannel::Rate &rate, const std::string &tab, std::size_t inst_seq, std::size_t gate_seq, CellInternalSignature::ComponentSubSignature &component){
						std::string for_what = "HHRate BaseRate "+itos(gate_seq)+" for Fixed channel "+itos(inst_seq);
						return 	DescribeRate_Thing( rate, tab, for_what, "r", component );
					};
					auto DescribeRate_Variable = [&DescribeRate_Thing](const IonChannel::Rate &rate, const std::string &tab, std::size_t inst_seq, std::size_t gate_seq, CellInternalSignature::ComponentSubSignature &component){
						std::string for_what = "HHRate BaseInf "+itos(gate_seq)+" for Fixed channel "+itos(inst_seq); // TODO more structured commenting
						return 	DescribeRate_Thing( rate, tab, for_what, "x", component );
					};
					auto DescribeRate_Tau = [&DescribeRate_Thing](const IonChannel::Rate &rate, const std::string &tab, std::size_t inst_seq, std::size_t gate_seq, CellInternalSignature::ComponentSubSignature &component){
						std::string for_what = "HHRate BaseTau "+itos(gate_seq)+" for Fixed channel "+itos(inst_seq); // TODO more structured commenting
						//printf("wac\n");
						return 	DescribeRate_Thing( rate, tab, for_what, "t", component );
					};
					
					auto DescribeRate_Q10 = [&](const Q10Settings &q10, const std::string & for_what, auto &pergate){
						char tmps[2000];
						if( q10.type == Q10Settings::FIXED ){
							pergate.Index_Q10 = AppendConstant( q10.q10, " Q10 Factor" );
							sprintf( tmps, "local_constants[%ld]", pergate.Index_Q10 );
							return std::string(tmps);
						}
						else if( q10.type == Q10Settings::FACTOR ){
							// NB equivalent to exp( ln(q10)*(temp - baseTemp)/10 )
							pergate.Index_Q10 = AppendConstant( q10.q10, " Q10 Factor" );
							pergate.Index_Q10_BaseTemp = AppendConstant( q10.experimentalTemp, " Q10 Base Temperature" );
							sprintf( tmps, "powf(local_constants[%ld], ( temperature - local_constants[%ld] ) / 10 )", pergate.Index_Q10, pergate.Index_Q10_BaseTemp );
							return std::string(tmps);
						}
						else return std::string("1");
					};
					
					// TODO actually let LEMS components receive rateScale perhaps? for model complenteness and nonlinearities ?
					// TODO also handle in each case, as q10 rate scaling arises
					// ccde +=   "		float rateScale = q10;"
					ccde +=   "		float rateScale = 1;\n";
					
					std::vector<std::string> factor_code_per_gates;
					for(std::size_t gate_seq = 0; gate_seq < chan.gates.contents.size(); gate_seq++){
						const auto &gate = chan.gates.contents[gate_seq];
						auto &pergate = distimpl.per_gate[gate_seq];
						const auto &for_that = for_what;
						std::string for_what = for_that+" channel "+itos(inst_seq);
						
						const std::string TauInf_suffix = Convert::Suffix(
							(( Scales<Time>::native ^ -1 ) * Scales<Time>::native ) // to unitless
						);
						
						auto Update_TauInf_Inline = [ &TauInf_suffix ]( auto Index_Q, const std::string &tab){
							std::string tauinf_code;
							char tmps[1000];
							
							tauinf_code +=   tab+"if(initial_state){\n";
							sprintf(tmps, "	local_stateNext[%ld] = inf;\n", Index_Q); tauinf_code += tab+tmps;
							tauinf_code +=   tab+"}else{\n";
							//sprintf(tmps, "	local_stateNext[%zd] = %s + dt * ( alpha * ( 1 - %s ) - beta * (%s) ) * q10 %s;\n", Index_Q, fana, fana, fana, Rate_suffix.c_str() ); ccde += tab+tmps;
							// sprintf(tmps, "	local_stateNext[%ld] = local_state[%ld] + dt * ( ( inf - local_state[%ld] ) / tau ) * q10 %s;\n", Index_Q, Index_Q, Index_Q, TauInf_suffix.c_str() ); tauinf_code += tab+tmps;
							
							// how many tau's does the timestep advance by?
							sprintf(tmps, "	float tau_factor = (( dt * q10)/ tau) %s;\n", TauInf_suffix.c_str() ); tauinf_code += tab+tmps;
							// blend current and inf values, bu a factor decaying from 1(now) to 0 (steady state)
							// tauinf_code +="	float blend_factor = 1/( 1 + tau_factor );\n"; // Bwd Euler (assuming dynamics to be fixed during dt)
							tauinf_code +="	float blend_factor = expf( -tau_factor );\n"; // Exp Euler (assuming dynamics to be fixed during dt). A.k.a. "cnexp"
							sprintf(tmps, "	local_stateNext[%ld] = (blend_factor) * local_state[%ld] + (1-blend_factor) * inf;\n", Index_Q, Index_Q ); tauinf_code += tab+tmps;
							// TODO in the case of nan's, prefer clipping to 0 or to 1 ?
							sprintf(tmps, "	if(!( local_stateNext[%ld] > (float)(1e-6) )) local_stateNext[%ld] = 1e-6;\n", Index_Q, Index_Q ); tauinf_code += tab+tmps;
							sprintf(tmps, "	if(!( local_stateNext[%ld] < (float)(1-1e-6) )) local_stateNext[%ld] = 1-1e-6;\n", Index_Q, Index_Q ); tauinf_code += tab+tmps;
							
							tauinf_code +=   "		}\n";
							
							return tauinf_code;
						};
						
						// Expose the gate variable, whatever its representation is supposed to be (single state, composition etc.)
						sprintf(tmps,"chan_gate_%zd_q", gate_seq);
						const std::string factor_name = tmps; const char *fana = factor_name.c_str();
						sprintf(tmps, "	float %s; \n", fana ); ccde += tmps;
						
						// Require amount of instances to be gathered from the specific sub-object containing the value
						Int instances = -1;
						
						if(gate.type == IonChannel::Gate::INSTANTANEOUS){
							auto instantaneous = gate.instantaneous;
							
							// Interface with gate rates
							ccde +=   "		{\n";
							ccde +=   		DescribeRate_Variable(instantaneous.steadyState, "\t\t", inst_seq, gate_seq, pergate.inf_component); 
							
							// Determine the gate variable
							sprintf(tmps, "		%s = x;\n", fana); ccde += tmps;
							ccde +=   "		}\n";
							
							instances = instantaneous.instances;
						}
						else if(
							gate.type == IonChannel::Gate::RATES
							|| gate.type == IonChannel::Gate::RATESTAU
							|| gate.type == IonChannel::Gate::RATESINF
							|| gate.type == IonChannel::Gate::RATESTAUINF
							|| gate.type == IonChannel::Gate::TAUINF
						){
							auto &gaga = gate.gaga;
							
							// in the origninal LEMS definitions, 'tau' or 'inf' terms are used if present, otherwise they are deduced from alpha, beta
							
							bool has_rates =
									gate.type == IonChannel::Gate::RATES
								||  gate.type == IonChannel::Gate::RATESTAU
								||  gate.type == IonChannel::Gate::RATESINF
								||  gate.type == IonChannel::Gate::RATESTAUINF
							; // otherwise it must be TauInf
							bool has_tau =
									gate.type == IonChannel::Gate::RATESTAU
								||  gate.type == IonChannel::Gate::RATESTAUINF
								||  gate.type == IonChannel::Gate::TAUINF
							; // otherwise it must have rates
							
							bool has_inf =
									gate.type == IonChannel::Gate::RATESINF
								||  gate.type == IonChannel::Gate::RATESTAUINF
								||  gate.type == IonChannel::Gate::TAUINF
							; // otherwise it must have rates
							
							float initial = NAN; // actually will be re-initialized at run time, to support LEMS components too
							
							if( has_tau ){
								initial       = DescribeRateThing::Value(comp_def.V0, gaga.steadyState);
							}
							else{
								double a_init = DescribeRateThing::Value(comp_def.V0, gaga.forwardRate);
								double b_init = DescribeRateThing::Value(comp_def.V0, gaga.reverseRate);
								
								initial = a_init / ( a_init + b_init );
							}
							
							
							// Allocate the gate state variable
							Int Index_Q = pergate.Index_Q = AppendStateVariable(initial, "Gatevar "+itos(gate_seq)+" for Fixed channel "+itos(inst_seq) );
							
							// Interface with gate rates
							
							// Determine the gate variable
							sprintf(tmps, "	%s = local_state[%ld]; \n", fana, Index_Q ); ccde += tmps;
							
							// add the internal dynamics code
							sprintf(tmps, "	// dynamics for channel %zd gate %zd \n", inst_seq, gate_seq); ccde += tmps;
							ccde +=   "	{\n";
							ccde +=   "		float q10 = "+DescribeRate_Q10(gaga.q10, factor_name, pergate)+";\n";
							
							if( has_rates ){
								ccde +=   "		float alpha;\n";
								ccde +=   "		{\n";
								ccde +=   		DescribeRate_Rate(gaga.forwardRate, "\t\t", inst_seq, gate_seq, pergate.alpha_component); 
								sprintf(tmps, "		alpha = r;\n"); ccde += tmps;
								ccde +=   "		}\n";
								
								ccde +=   "		float beta;\n";
								ccde +=   "		{\n";
								ccde +=   		DescribeRate_Rate(gaga.reverseRate, "\t\t", inst_seq, gate_seq, pergate.beta_component);
								sprintf(tmps, "		beta = r;\n"); ccde += tmps;
								ccde +=   "		}\n";
								if(config.debug){
									// ccde += "	printf(\"albe %g %g\\n\", alpha, beta);\n";
								}
							}
							
							ccde +=   "		float tau;\n";
							if( has_tau ){
								ccde +=   "		{\n";
								ccde +=   		DescribeRate_Tau(gaga.timeCourse, "\t\t", inst_seq, gate_seq, pergate.tau_component); 
								sprintf(tmps, "		tau = t;\n"); ccde += tmps;
								ccde +=   "		}\n";
							}
							else{
								ccde +=   "		tau = 1 / ( alpha + beta );\n";
							}
							
							ccde +=   "		float inf;\n";
							if( has_inf ){
								ccde +=   "		{\n";
								ccde +=   		DescribeRate_Variable(gaga.steadyState, "\t\t", inst_seq, gate_seq, pergate.inf_component); 
								sprintf(tmps, "		inf = x;\n"); ccde += tmps;
								ccde +=   "		}\n";
							}
							else{
								ccde +=   "		inf = alpha / ( alpha + beta );\n";
							}
							
							
							ccde += Update_TauInf_Inline( Index_Q, "\t\t");
							
							// end internal dynamics
							ccde +=   "	}\n";
							
							instances = gaga.instances;
						}
						else if( gate.type == IonChannel::Gate::FRACTIONAL ){
							const auto &fga = chan.fractional_gates.at( gate.fractional );
							
							pergate.Index_Q = -1; // composite
							
							// Determine the gate variable
							ccde += "	" + std::string(fana) + " = 0;\n"; // sum of subgate vars
							
							// Also define q10, common for all 
							sprintf(tmps, "	// dynamics for %s \n", for_what.c_str() ); ccde += tmps;
							ccde +=   "	{\n";
							ccde +=   "		float q10 = "+DescribeRate_Q10(fga.q10, factor_name, pergate)+";\n";
							
							for( size_t sga_seq = 0; sga_seq < fga.subgates.size() ; sga_seq++ ){
								const auto &sga = fga.subgates.at(sga_seq);
								
								CellInternalSignature::IonChannelDistImplementation::SubGate persub;
								
								const std::string &for_that = for_what;
								std::string for_what = for_that+" subgate "+itos(sga_seq);
								
								// Allocate the subgate state variable
								float initial = DescribeRateThing::Value(comp_def.V0, sga.steadyState);
								Int Index_SubQ = persub.Index_Q = AppendStateVariable( initial, for_what + " Variable" );
								
								// and its contribution factor
								std::size_t Index_SubQFactor = AppendStateVariable( sga.fraction_of_conductivity, for_what + " Effective Fraction" );
								
								// contribute to gate variable
								sprintf(tmps, "	%s += local_state[%ld] * local_constants[%zd]; \n", fana, Index_SubQ, Index_SubQFactor ); ccde += tmps;
								
								// Interface with subgate rates
								
								// add the internal dynamics code
								sprintf(tmps, "	// dynamics for %s \n", for_what.c_str()); ccde += tmps;
								ccde +=   "	{\n";
								// subgate Q10 is defined in the schema but not used, somehow ???
								
								ccde +=   "		float tau;\n";
								ccde +=   "		{\n";
								ccde +=   		DescribeRate_Tau(sga.timeCourse, "\t\t", inst_seq, gate_seq, persub.tau_component); 
								sprintf(tmps, "		tau = t;\n"); ccde += tmps;
								ccde +=   "		}\n";
								
								ccde +=   "		float inf;\n";
								ccde +=   "		{\n";
								ccde +=   		DescribeRate_Variable(sga.steadyState, "\t\t", inst_seq, gate_seq, persub.inf_component); 
								sprintf(tmps, "		inf = x;\n"); ccde += tmps;
								ccde +=   "		}\n";
								
								
								ccde += Update_TauInf_Inline( Index_SubQ, "\t\t");
							
								pergate.subgates.push_back(persub);
								
								// end internal dynamics
								ccde +=   "	}\n";
								
								
							}
							ccde +=   "	}\n";
							
							instances = fga.instances;
						}
						else if( gate.type == IonChannel::Gate::KINETIC ){
							const auto &ks = chan.kinetic_gates.at( gate.kinetic );
							
							const size_t states = ks.state_names.size();
							
							pergate.Index_Q = -1; // kinetic
							
							// Also define q10, common for all 
							sprintf(tmps, "	// dynamics for %s \n", for_what.c_str()); ccde += tmps;
							ccde +=   "	{\n";
							ccde +=   "		float q10 = "+DescribeRate_Q10(ks.q10, factor_name, pergate)+";\n";
							
							// declare state variables
							for( size_t state_seq = 0; state_seq < states; state_seq++ ){
								CellInternalSignature::IonChannelDistImplementation::SubGate persub;
								
								const std::string &for_that = for_what;
								std::string for_what = for_that+" state "+itos(state_seq);
								persub.Index_Q = AppendStateVariable( NAN, for_what + " Variable" );
								
								pergate.subgates.push_back(persub);
							}
							// gather transition rates
							struct FromToInfo{
								Int tran_seq;
								Int from;
								Int to;
							};
							std::vector< std::string > transition_names( ks.transitions.size() );
							std::vector< std::vector< FromToInfo > > trans_from(states), trans_to(states); // per state
							
							for( Int tran_seq = 0; tran_seq < (Int)ks.transitions.size(); tran_seq++ ){
								auto &transition = ks.transitions[tran_seq];
								
								CellInternalSignature::IonChannelDistImplementation::SubGate pertran;
								Int from, to;
								std::string ratecode;
								
								sprintf(tmps, "	// dynamics for transition %ld \n", tran_seq); ccde += tmps;
								
								ratecode += "	float alpha, beta;\n";
								ratecode += "	{\n";
								if( transition.type == IonChannel::GateKS::Transition::FORWARD_REVERSE ){
									auto &forrev = transition.forrev;
									from = forrev.from; to = forrev.to;
									
									ratecode +=   "		{\n";
									ratecode +=   		DescribeRate_Rate(forrev.forwardRate, "\t\t", inst_seq, gate_seq, pergate.alpha_component); 
									sprintf(tmps, "		alpha = r;\n"); ratecode += tmps;
									ratecode +=   "		}\n";
									
									ratecode +=   "		{\n";
									ratecode +=   		DescribeRate_Rate(forrev.reverseRate, "\t\t", inst_seq, gate_seq, pergate.beta_component); 
									sprintf(tmps, "		beta = r;\n"); ratecode += tmps;
									ratecode +=   "		}\n";
									
									
								}
								else if( transition.type == IonChannel::GateKS::Transition::TAU_INF ){
									auto &tauinf = transition.tauinf;
									from = tauinf.from; to = tauinf.to;
									
									ratecode +=   "		float tau;\n";
									ratecode +=   "		{\n";
									ratecode +=   		DescribeRate_Tau(tauinf.timeCourse, "\t\t", inst_seq, gate_seq, pertran.tau_component); 
									sprintf(tmps, "		tau = t;\n"); ratecode += tmps;
									ratecode +=   "		}\n";
									
									ratecode +=   "		float inf;\n";
									ratecode +=   "		{\n";
									ratecode +=   		DescribeRate_Variable(tauinf.steadyState, "\t\t", inst_seq, gate_seq, pertran.inf_component); 
									sprintf(tmps, "		inf = x;\n"); ratecode += tmps;
									ratecode +=   "		}\n";
									
									ratecode +=   "		alpha = inf / tau;\n";
									ratecode +=   "		beta  = ( 1 - inf ) / tau;\n";
									
								}
								else{
									printf("ks implementation transition type\n");
									return false;	
								}
								ratecode += "	}\n";
								
								auto tranname = transition_names[tran_seq] = "transition_from_"+itos(from)+"_to_"+itos(to);
								
								ccde +=   "	float "+tranname+"_for"+" = NAN;\n";
								ccde +=   "	float "+tranname+"_rev"+" = NAN;\n";
								ccde +=   "	{\n";
								
								
								ccde +=   ratecode;
								
								ccde +=   "	"+tranname+"_for"+" = alpha;\n";
								ccde +=   "	"+tranname+"_rev"+" = beta;\n";
								ccde +=   "	}\n";
								
								trans_from[from].push_back( { tran_seq, from, to} );
								trans_to  [to  ].push_back( { tran_seq, from, to} );
								
								// in this case, fluxes are bidirectional (at least in structure)
								trans_to  [from].push_back( { tran_seq, from, to} );
								trans_from[to  ].push_back( { tran_seq, from, to} );
								
								pergate.transitions.push_back(pertran);
								
							}
							// update
							ccde +=   "	if(initial_state){\n";
							ccde +=   "		// XXX no initial constants specified :(\n"
										"		// init to all to first state\n";
							// XXX NeuroML and LEMS do not specify initialization of kinetic schemes
							// LEMS initializes to all quantity to state 0 :
							//  https://github.com/LEMS/jLEMS/blob/af575ea517f446bb4fd7389c427eb418e1a46d86/src/main/java/org/lemsml/jlems/core/run/KScheme.java#L70
							// perhaps try solving the matrix (assuming all rates are constant, and matrix is invertible) LATER
							ccde +=   "		//FIXME\n";
							
							for( size_t state_seq = 0; state_seq < states; state_seq++ ){
								auto Index_Q = pergate.subgates.at(state_seq).Index_Q;
								sprintf(tmps, "			local_stateNext[%ld] = %d;\n", Index_Q, (state_seq == 0) ? 1 : 0 ); ccde += tmps;
							}
							ccde +=   "	}else{\n";
							
							for( size_t state_seq = 0; state_seq < states; state_seq++ ){
								
								// gather rates flowing to this state variable 
								
								sprintf(tmps, "			float flux_offdiag_%zd = 0", state_seq ); ccde += tmps;
								for( auto tranto   : trans_to  .at(state_seq) ){
									
									// now, the flux toward this state
									// could either be by the 'forward' or the 'reverse' part of the transition.
									
									// If 'from' is not state_seq, then it's the forward rate, toward state_seq
									const char *sDirection = "for";
									auto actual_from = tranto.from;
									// otherwise, 'to' is flowing to 'from', with the reverse rate
									if( actual_from == (Int)state_seq ){
										// the reverse flux of the transition
										sDirection = "rev";
										actual_from = tranto.to;
									}
									
									sprintf(tmps, "	+ ( %s_%s * local_state[%ld] )",
										transition_names.at(tranto.tran_seq).c_str(),
										sDirection,
										pergate.subgates.at(actual_from).Index_Q
									); ccde += tmps;
								}
								ccde += ";\n";
								
								// gather rates flowing from this state variable
								sprintf(tmps, "			float rate_diag_%zd = 0", state_seq ); ccde += tmps;
								for( auto tranfrom   : trans_from  .at(state_seq) ){
									
									// determine direction as with off-diagonal fluxes
									const char *sDirection = "for";
									// otherwise, 'to' is flowing to 'from', with the reverse rate
									if( tranfrom.from != (Int)state_seq ){
										// the reverse flux of the transition
										sDirection = "rev";
									}
									
									sprintf(tmps, "	+ %s_%s",
										transition_names.at(tranfrom.tran_seq).c_str(),
										sDirection
									); ccde += tmps;
								}
								ccde += ";\n";
								
							}
							for( size_t state_seq = 0; state_seq < states; state_seq++ ){
								auto Index_Q = pergate.subgates.at(state_seq).Index_Q;
								sprintf(tmps, "			local_stateNext[%ld] = local_state[%ld] + dt * ( flux_offdiag_%zd - ( rate_diag_%zd * local_state[%ld] ) ", Index_Q, Index_Q, state_seq, state_seq, Index_Q ); ccde += tmps;
								// for( auto tranto   : trans_to  .at(state_seq) ){
								// 	sprintf(tmps, "	+ ( %s_rev * local_state[%ld] ) ", transition_names.at(tranto.tran_seq).c_str(), pergate.subgates.at(tranto.to).Index_Q ); ccde += tmps;
								// }
								// sprintf(tmps, "	- ( local_stateNext[%ld] * ( 0", Index_Q ); ccde += tmps;
								// for( auto tranfrom : trans_from.at(state_seq) ){
								// 	sprintf(tmps, "	+ %s_for ", transition_names.at(tranfrom.tran_seq).c_str() ); ccde += tmps;
								// }
								// ccde += " ) )"; // end transitions leaking from this state
								sprintf(tmps, " ) * q10 %s;\n", TauInf_suffix.c_str() ); ccde += tmps;
								
								sprintf(tmps, "			// add some sanity clipping\n" ); ccde += tmps;
								sprintf(tmps, "			if( local_stateNext[%ld] > 1 ) local_stateNext[%ld] = 1;\n", Index_Q, Index_Q ); ccde += tmps;
								sprintf(tmps, "			if( local_stateNext[%ld] < 0 ) local_stateNext[%ld] = 0;\n", Index_Q, Index_Q ); ccde += tmps;
								
							}
							sprintf(tmps, "			// finally, preserve a total of 1, divergence goes to first state as in NEURON\n" ); ccde += tmps;
							{
							auto Index_Q = pergate.subgates.at(0).Index_Q;
							sprintf(tmps, "			local_stateNext[%ld] = 1", Index_Q ); ccde += tmps;
							for( size_t state_seq = 0; state_seq < states; state_seq++ ){
								
								if( state_seq == 0 ) continue;
								auto Index_Q = pergate.subgates.at(state_seq).Index_Q;
								sprintf(tmps, " - local_stateNext[%ld]", Index_Q ); ccde += tmps;
							}
							sprintf(tmps, ";\n" ); ccde += tmps;
							
							}
							
							ccde +=   "		}\n";
							
							// Determine the gate variable
							ccde += "	" + std::string(fana) + " = 0";
							for( auto open : ks.open_states ){
								ccde += " + local_state[" + itos( pergate.subgates.at(open).Index_Q ) + "]";
							}
							ccde += ";\n"; // sum of subgate vars
							
							ccde +=   "	}\n";
							
							instances = ks.instances;
							
							// done !
						}
						else if( gate.type == IonChannel::Gate::COMPONENT ){
							ccde +=   "		{\n";
							// ccde += "		// LEMS gate component\n";
							ccde += DescribeLemsInline.SingleInstance( gate.component, "\t\t", for_what, pergate.inf_component, config.debug );
							sprintf(tmps, "		%s = Lems_exposure_fcond;\n", fana); ccde += tmps;
							ccde +=   "		}\n";
						}
						else{
							printf("internal error: odd ion channel gates not supported yet\n");
							return false;
						}
						
						// and use the gate variable
						if( gate.type == IonChannel::Gate::COMPONENT ){
							factor_code_per_gates.push_back( std::string("* ") + fana ); // fcond is fana
						}
						else{
							if( instances < 0 ){
								printf("internal error: instance count\n");
								return false;
							}
							std::string factor_code;
							for(int i = 0; i < instances; i++){
								factor_code += " * "; factor_code += fana;
							}
							factor_code_per_gates.push_back(factor_code); // fcond is (fana ** instances)
						}
						
					}
					
					std::string factor_string;
					for(size_t i = 0; i < factor_code_per_gates.size(); i++){
						factor_string += factor_code_per_gates[i];
					}
					
					// add conductance scaling
					ccde   += "	float conductance_scaling = 1;\n";
					if( chan.conductance_scaling.type == IonChannel::ConductanceScaling::NONE ){
						// it's ok, do nothing
					}
					else if( chan.conductance_scaling.type == IonChannel::ConductanceScaling::Q10 ){
						sprintf(tmps, "	conductance_scaling = %s;\n", DescribeRate_Q10( chan.conductance_scaling.q10, for_what + " Scaling Factor", distimpl.conductance_scaling ).c_str() ); ccde += tmps;
					}
					else if( chan.conductance_scaling.type == IonChannel::ConductanceScaling::COMPONENT ){
						ccde +=   "	{\n";
						// ccde += "	// LEMS gate conductance scaling\n";
						ccde += DescribeLemsInline.SingleInstance( chan.conductance_scaling.component, "\t", for_what, distimpl.conductance_scaling.scaling_component, config.debug );
						sprintf(tmps, "	conductance_scaling = Lems_exposure_factor;\n"); ccde += tmps;
						ccde +=   "	}\n";
					}
					else{
						printf("unknown conductance scaling type\n");
						return false;
					}
					
					sprintf(tmps, "		ChannelOpenFraction = conductance_scaling %s;\n", factor_string.c_str()); ccde += tmps;
				}
				// NB: Assume that ∂/∂V ChannelOpenFraction = 0 and take the associated risk of instability. Improve LATER when the need arises, because re-evaluating the whole ion channel for a delta is too much effort atm.
				
				// now determine channel distribution current (also the point to add LEMS channel distributions LATER)
				ccde   += "	float I_chan = NAN;\n";
				ccde   += "	float channel_conductance = NAN;\n";
				// a second pass through channel distributions, to determine what will be done with channel fopen exposure (GHK1 does something... different that could be seen as variable Gbase)
				if( provides_current ){
					if( is_population ){
						std::size_t Index_Channel_Number = AppendConstant( inst.number, for_what + " Population Count" );
						sprintf(tmps, "		float Population_Count = local_constants[%zd]; // conductivity\n", Index_Channel_Number); ccde += tmps;
						
						const std::string gPop_suffix = Convert::Suffix(
							Scales<Conductance>::native
							.to(Scales<Conductance>::native)
						);
						// FIXME ChannelConductance is only set if the ion channel is a omponent ... make sure that's right !
						sprintf(tmps, "		float gTotal = (Population_Count * ChannelConductance * ChannelOpenFraction)%s; //total conductance\n", gPop_suffix.c_str()); ccde += tmps;
						
						const std::string Ichan_suffix = Convert::Suffix(
							(Scales<Voltage>::native * Scales<Conductance>::native)
							.to(Scales<Current>::native)
						);
						sprintf(tmps, "		I_chan = ( gTotal * (Erev - Vcomp) )%s; //total current\n", Ichan_suffix.c_str()); ccde += tmps;
						ccde +=  "	`channel_conductance = gTotal;\n";
					}
					else{
						printf("internal error: uses what sort of current?\n");
						return false;
					}
				}
				else if( provides_density ){
					
					ccde   += "	float iDensity = NAN;\n";
					ccde   += "	float gDensity = NAN;\n";
					
					if( inst.type == ChannelDistribution::GHK ){
						double R = 8.3144621; // ( J / ( K * mol ) )
						double zCa = 2; // any ion LATER
						double F = 96485.3; // ( Cb / mol )
						
						std::size_t Index_Channel_Permeability = AppendConstant( inst.permeability,  for_what + " Permeability " );
						sprintf(tmps, "	float permeability = local_constants[%zd];\n", Index_Channel_Permeability); ccde += tmps;
						
						auto SI_To_InvVolt_suffix = Convert::Suffix( Scales<Dimensionless>::native.to( (Scales<Voltage>::native)^(-1) ) ); // dimensionless Volt units to Erev
					
						sprintf(tmps, "	float K = ( ( %.17g * %.17g) / (%.17g * temperature) )%s;\n", zCa, F, R, SI_To_InvVolt_suffix.c_str()); ccde += tmps;
						
						ccde   += "	float expKv = expf( -1 * K * Vcomp );\n";
						
						const std::string iDensity_suffix = Convert::Suffix(
							( Scales<Permeability>::native * Scales<Concentration>::native ) 
							.to( Scales<Current>::native / (microns^2) )
						);
						const std::string gDensity_suffix = Convert::Suffix(
							( (Scales<Current>::native / (microns^2)) / Scales<Voltage>::native ) 
							.to( Scales<Conductance>::native / (microns^2) )
						);// through explicit differentiation
						auto EvalJ = [&tmps, &iDensity_suffix, zCa, F ]( std::string &ccde ){
							sprintf(tmps, "		iDensity = (-1 * permeability * ChannelOpenFraction * %.17g * %.17g * K * Vcomp * ( Ca_concentration - (Ca_concentration_extra * expKv) ) / (1 - expKv))%s;\n", zCa, F, iDensity_suffix.c_str()); ccde += tmps;
						};
						
						ccde   += "	if( Ca_concentration_extra > 0 ){\n";
						
						// TODO make the delta V some sort of named constant
						ccde += "	const float delta_volt = ((float)(1e-6 "+Convert::Suffix(Scales<Dimensionless>::native.to( Scales<Voltage>::native ))+"));\n"; 
						ccde += "	float iDensity_at_plus_1uv = NAN;\n";
						ccde += "	{\n";
						ccde += "	float V_tmp = Vcomp; float Vcomp = V_tmp + delta_volt;\n";
						EvalJ( ccde );
						ccde += "	iDensity_at_plus_1uv = iDensity;\n";
						ccde += "	}\n";
						EvalJ( ccde );
						ccde += "	\n";
						
						// differentiate iDensity by sampling a microvolt apart, differentiate properly LATER if the need arises
						ccde   += "		gDensity = ( ( iDensity_at_plus_1uv - iDensity ) / (delta_volt) )"+gDensity_suffix+";\n"; // conductivity * ( volts / volts ) to ( current / sq.um )
						
						ccde   += "	}else{\n";
						ccde   += "		iDensity = 0;\n";
						ccde   += "		gDensity = 0;\n";
						ccde   += "	}\n";
					}
					else if( uses_conductivity ){
						// it is a conductivity-based thing
						const std::string iDensity_suffix = Convert::Suffix(
							( Scales<Voltage>::native * Scales<Conductivity>::native ) 
							.to( Scales<Current>::native / (microns^2) )
						);
						const std::string gDensity_suffix = Convert::Suffix(
							( Scales<Conductivity>::native ) 
							.to( Scales<Conductance>::native / (microns^2) )
						);
						
						if( 
							uses_fixed_conductivity 
						){
							std::size_t Index_Channel_Gbase = AppendConstant( inst.conductivity, for_what + " Total Base Conductivity" );
							sprintf(tmps, "		float Gbase = local_constants[%zd]; // conductivity\n", Index_Channel_Gbase); ccde += tmps;
						}
						else{
							printf(" internal error: ion channel distribution with conductivity and no Gbase\n");
							return false;
						}
						
						// apply the scaling factor
						sprintf(tmps, "		float Gscaled  = Gbase * ChannelOpenFraction;\n"); ccde += tmps;
						
						// does it modify Gbase somehow?
						// NB: in the case of GHK2, popen is voltage actually
						if( inst.type == ChannelDistribution::GHK2 ){
							const std::string UnitToVolt_suffix = Convert::Suffix(
								Scales<Dimensionless>::native 
								.to( Scales<Voltage>::native )
							);
								// modify Gbase
							ccde   += " float tmp = ( 25 * temperature ) / ( 293.15 * 2 ); // unitless kelvins\n";
							
							auto EvalPopen = [&tmps, &UnitToVolt_suffix ]( std::string &ccde ){
								sprintf(tmps, "	float V = Vcomp * ( 1000 / (1%s) ); // unitless millivolts\n", UnitToVolt_suffix.c_str() ); ccde += tmps;
								ccde   += " float pOpen = NAN;\n";
								ccde   += "	if( Vcomp == 0 ){\n"; // TODO perhaps add tolerance around some epsilon
								ccde   += "		pOpen = tmp * ( 1 - ( Ca_concentration / Ca_concentration_extra ) ) * (1e-3 "+UnitToVolt_suffix+");\n";
								ccde   += "	}else{\n";
								ccde   += "		pOpen = tmp * ( 1 - ( ( Ca_concentration / Ca_concentration_extra ) * expf( V / tmp ) ) ) * ( ( V / tmp ) / ( expf( V / tmp ) - 1) ) * (1e-3"+UnitToVolt_suffix+");\n";
								ccde   += "	}\n";
								
								ccde   += "	if( Ca_concentration_extra == 0 ){\n";
								ccde   += "		pOpen = 0;\n";
								ccde   += "	}\n";
							};
							
							// TODO make the delta V some sort of named constant
							ccde += "	const float delta_volt = ((float)(1e-6 "+Convert::Suffix(Scales<Dimensionless>::native.to( Scales<Voltage>::native ))+"));\n"; 
							ccde += "	float pOpen_at_plus_1uv = NAN;\n";
							ccde += "	{\n";
							ccde += "	float V_tmp = Vcomp; float Vcomp = V_tmp + delta_volt;\n";
							EvalPopen( ccde );
							ccde += "	pOpen_at_plus_1uv = pOpen;\n";
							ccde += "	}\n";
							EvalPopen( ccde );
							ccde += "	\n";
							
							// differentiate pOpen by sampling a microvolt apart, differentiate properly LATER if the need arises
							ccde   += "	iDensity = (Gscaled  * pOpen)"+iDensity_suffix+";\n";
							ccde   += "	gDensity = Gscaled * ( ( pOpen_at_plus_1uv - pOpen ) / (delta_volt) )"+gDensity_suffix+";\n"; // conductivity * ( volts / volts ) to ( current / sq.um )
							// ccde   += "	printf(\"popen %f\\n\", pOpen);\n";
							// assuming Gscaled is constant wrt Vcomp
						}
						else{
							sprintf(tmps, "	iDensity = Gscaled * (Erev - Vcomp)%s;\n", iDensity_suffix.c_str()); ccde += tmps;
							ccde   += "	gDensity = (Gscaled)"+gDensity_suffix+";\n";
						}
						
					}
					else{
						printf("internal error: ion channel uses what sort of current density?\n");
						return false;
					}
					
					// finally, total compartment current
					const std::string Ichan_suffix = Convert::Suffix(
						(	
							( Scales<Current>::native / (microns^2) ) 
							* (microns^2)
						)
						.to( Scales<Current>::native )
					);
					const std::string Gchan_suffix = Convert::Suffix(
						(	
							( Scales<Conductance>::native / (microns^2) ) 
							* (microns^2)
						)
						.to( Scales<Conductance>::native )
					);
					sprintf(tmps, "	I_chan = (iDensity * Acomp)%s;\n", Ichan_suffix.c_str()); ccde += tmps;
					sprintf(tmps, "	channel_conductance = (gDensity * Acomp)%s;\n", Gchan_suffix.c_str()); ccde += tmps;
				}
				else{
					// component populations MUCH LATER ?
					printf("internal error: ion channel provides what, if not current or density?\n");
					return false;
				}
				
				// contribute to compartment current influx
				sprintf(tmps, "		I_channels_total += I_chan;\n"); ccde += tmps;
				
				// sprintf(tmps, "		assert(isfinite(channel_conductance));\n"); ccde += tmps; should be done some time...
				sprintf(tmps, "		G_channels_total += channel_conductance;\n"); ccde += tmps;
				
				// also contribute to ion channels
				Int species = chan.species; 
				if( comp_def.ions.count(species) > 0 ){
					sprintf(tmps, "		I_ion_%ld += I_chan;\n", species); ccde += tmps;
				}
				// TODO g_ion as well, to improve integration stability
				
				ccde   += "\n";
				
				ccde   += "	}\n";
				ccde   += "\n";
				
				//printf("(done channel inst %zd)\n", inst_seq);
			}
			
			// add synapses
			ccde += "	// Current from synapses\n";
			sprintf(tmps, "	float I_synapses_total = 0;\n"); ccde += tmps;
			sprintf(tmps, "	float G_synapses_total = 0;\n"); ccde += tmps; // FIXME
			
			for(Int id_id : comp_def.synaptic_component_types.toArray()){
				if( !ImplementSynapseType( AppendSingle, AppendMulti, DescribeLemsInline, random_call_counter, for_what + " Synapse type "+std::to_string(id_id), tab, id_id, comp_impl.synapse, ccde ) ) return false;
			}
			
			// add inputs
			ccde += "	// Current from inputs\n";
			sprintf(tmps, "	float I_input_total = 0;\n"); ccde += tmps;
			sprintf(tmps, "	float G_input_total = 0;\n"); ccde += tmps;
			
				
			// generate tables and code for each input type
			for(Int id_id : comp_def.input_types.toArray()){
				if( !ImplementInputSource( AppendSingle, AppendMulti, DescribeLemsInline, random_call_counter, for_what + " Input type "+std::to_string(id_id), tab, id_id, comp_impl.input, ccde ) ) return false;
			}
			
			// all ion contributions have been considered; now integrate the ion dynamics
			for( auto keyval : comp_def.ions){
				Int species_seq = keyval.first;
				const auto &instance = keyval.second;
				
				// now work with ion pool :D
				const CellInternalSignature::IonSpeciesDistImplementation &distimpl = comp_impl.concentration.at(species_seq);
				
				const auto &for_that = for_what;
				std::string for_what = for_that+" Ion "+itos(species_seq)+" pool";
				
				const auto &conc_model = conc_models.get(instance.conc_model_seq);
				
				const std::string tab = "\t";
				std::string ionpool_code;
				char tmps[2000];
				
				
				sprintf(tmps, "	// Dynamics for ion %ld \n", species_seq); ionpool_code += tmps;
				ionpool_code += tab + "{\n";
				ionpool_code += tab;
				
				// TODO deduplicate also in dynamics section
				ionpool_code += ExposeRequirements_ConcModel( species_seq, distimpl, tab );
				
				// no compile time dimensions for moles and coulombs, but they cancel each other with the Faraday constant so it all checks out
				
				
				if(conc_model.type == ConcentrationModel::COMPONENT){
					
					ionpool_code +=  tab+"// LEMS component\n";
					const auto &comptype = model.component_types.get(conc_model.component.id_seq);
					
					std::string lemscode = DescribeLems::Assigned(comptype, model.dimensions, distimpl.component, &AppendSingle, for_what, tab, random_call_counter, config.debug);
					ionpool_code += lemscode;
					
					// numerical integration code here
					std::string lemsupdate = DescribeLems::Update(comptype, model.dimensions, distimpl.component, &AppendSingle, for_what, tab, random_call_counter, config.debug );
					ionpool_code += lemsupdate;
					
					ionpool_code += DescribeLems::Exposures(comptype, for_what, tab, config.debug);
				}
				else{
					// it is a built-in type
					sprintf(tmps, " float iCa = I_ion_%ld; //total ion current\n", keyval.first); ccde += tmps;
					
					ionpool_code += tab+"float ion_charge = 2;\n"; //TODO
					
					ionpool_code +=  tab+"float influx_rate = NAN;\n";
					if(conc_model.type == ConcentrationModel::LEAKY){
						const std::string CurrentToConcRate_suffix = Convert::Suffix(
							( Scales<Current>::native / (microns^3) )
							.to( Scales<Concentration>::native / Scales<Time>::native )
						);
						
						sprintf(tmps, "float Faraday = %.17g;\n", 96485.3); ionpool_code += tab+tmps;
						
						sprintf(tmps, "float shellThickness = local_constants[%zd];\n", distimpl.Index_Shellthickness_Or_RhoFactor); ionpool_code += tab+tmps;
						// TODO check dimensions & units here !
						sprintf(tmps, "float effectiveRadius = sqrt(Acomp / (float)(4 * M_PI));\n" ); ionpool_code += tab+tmps;
						sprintf(tmps, "float innerRadius = effectiveRadius - shellThickness;\n" ); ionpool_code += tab+tmps;
						sprintf(tmps, "float shellVolume = ((effectiveRadius * effectiveRadius * effectiveRadius) * (float)(4 * M_PI / 3)) - ((innerRadius * innerRadius * innerRadius) * (float)(4 * M_PI / 3));\n"); ionpool_code += tab+tmps;
						sprintf(tmps, "influx_rate = ( iCa / (ion_charge * Faraday * shellVolume) )%s;\n", CurrentToConcRate_suffix.c_str()); ionpool_code += tab+tmps;
						// if( config.debug ){
						// 	ionpool_code += "		printf(\"effectiveRadius %e \\ninnerRadius %e\\nshellVolume %e\\n\", effectiveRadius, innerRadius, shellVolume);\n";
						// }
					}
					else if(conc_model.type == ConcentrationModel::FIXED_FACTOR){
						const std::string CurrentToConcRate_suffix = Convert::Suffix(
							( (Scales<Current>::native / (microns^2)) * Scales<RhoFactor>::native  )
							.to( Scales<Concentration>::native / Scales<Time>::native )
						);
						
						sprintf(tmps, "influx_rate = ( (iCa / Acomp) * local_constants[%zd] )%s;\n", distimpl.Index_Shellthickness_Or_RhoFactor, CurrentToConcRate_suffix.c_str()); ionpool_code += tab+tmps;
					}
					else{
						assert(false); // LATER something better for internal errors
					}
					
					ionpool_code += tab+"if(initial_state){\n";
					ionpool_code += tab+"	// initialize\n";
					sprintf(tmps, "		local_stateNext[%zd] = local_state[%zd];\n", distimpl.Index_Intra,  distimpl.Index_Intra); ionpool_code += tab+tmps;
					sprintf(tmps, "		local_stateNext[%zd] = local_state[%zd];\n", distimpl.Index_Extra,  distimpl.Index_Extra); ionpool_code += tab+tmps;
					ionpool_code += "	}else{\n";
					
					const std::string ConcToConcRate_suffix = Convert::Suffix(
						( Scales<Concentration>::native / Scales<Time>::native )
						.to( Scales<Concentration>::native / Scales<Time>::native )
					);
					
					sprintf(tmps, "		float leak_rate = ( ( local_state[%zd] - local_constants[%zd] ) / local_constants[%zd] )%s;\n", distimpl.Index_Intra, distimpl.Index_RestConc, distimpl.Index_DecayTau, ConcToConcRate_suffix.c_str()); ionpool_code += tab+tmps;
					// if(config.debug){
					// ionpool_code += "\t\tprintf(\"I %e In %e Leak %e\\n\", iCa, influx_rate, leak_rate);\n";
					// //sprintf(tmps, "\t\tprintf(\"%%g\\t%%g\\t%%g\\t%%g\\n\", I_axial, I_channels_total, I_input_total, I_synapses_total);\n"); ionpool_code += tmps;
					// }
					sprintf(tmps, "		local_stateNext[%zd] = local_state[%zd] + ( dt * ( influx_rate - leak_rate ) );\n", distimpl.Index_Intra, distimpl.Index_Intra); ionpool_code += tab+tmps;
					sprintf(tmps, "		if( local_stateNext[%zd] < 0 ) local_stateNext[%zd] = 0;\n", distimpl.Index_Intra, distimpl.Index_Intra); ionpool_code += tab+tmps;
					// no changes for extra
					sprintf(tmps, "		local_stateNext[%zd] = local_state[%zd];\n", distimpl.Index_Extra,  distimpl.Index_Extra); ionpool_code += tab+tmps;
					ionpool_code += tab+"}\n";
					
				}
				
				ionpool_code += tab + "}\n";
				ccde += ionpool_code;
				
			}
			
			// finally, integrate currents into voltage
			// TODO the cable equation handling should be broken off since it involves other compartments, therefore it is not really internal (important for compartment as work itam)
			sprintf(tmps, "	I_internal = I_channels_total + I_input_total + I_synapses_total;\n"); ccde += tmps;
			sprintf(tmps, "	G_internal = G_channels_total + G_input_total + G_synapses_total;\n"); ccde += tmps;
			
			// TODO compute axial conductance regardless of cable solver
			
			// store diagonal of cable equation matrix
			std::string InvRC_suffix = Convert::Suffix( ( (Scales<Conductance>::native / Scales<Capacitance>::native) ).to( Scales<Frequency>::native ) );
			if( cell_cable_solver == SimulatorConfig::CABLE_BWD_EULER ){
				sprintf(tmps, "		PerComp_InvRC_Diagonal[comp] = 1*(G_internal / C[comp] )%s + PerComp_InvRC_Axial[comp];\n", InvRC_suffix.c_str()); ccde += tmps;
			}
			ccde += "\n";
			
			
			ccde += "	if(initial_state){\n";
			ccde += "		// initialize\n";
			sprintf(tmps, "		V_next[comp] = V[comp];\n"); ccde += tmps;
			ccde += "	}else{\n";
			if(config.debug){
			ccde += "		printf(\"axial\\tchannel\\tinput\\tsyn\\n\");\n";
			if( uses_Iaxial ){
				sprintf(tmps, "		printf(\"iax %%g\\t%%g\\t%%g\\t%%g\\n\", I_axial, I_channels_total, I_input_total, I_synapses_total);\n"); ccde += tmps;
			}
			else{
				sprintf(tmps, "		printf(\"noiax %%g\\t%%g\\t%%g\\n\", I_channels_total, I_input_total, I_synapses_total);\n"); ccde += tmps;
				sprintf(tmps, "		printf(\"nogax %%g -- %%g\\t%%g\\t%%g\\n\", PerComp_InvRC_Diagonal[comp], G_channels_total, G_input_total, G_synapses_total);\n"); ccde += tmps;
			}
			}
			const std::string Vnext_suffix = Convert::Suffix(
				(Scales<Time>::native * Scales<Current>::native / Scales<Capacitance>::native)
				.to(Scales<Voltage>::native)
			);
			
			// axial currents are integrated otherwise, in backward(ish) cable equation solvers
			if( cell_cable_solver == SimulatorConfig::CABLE_FWD_EULER ){
				sprintf(tmps, "		V_next[comp] = V[comp] + ( dt * ( I_internal + I_axial ) / C[comp] )%s;\n", Vnext_suffix.c_str()); ccde += tmps;
			}
			else{
				// sprintf(tmps, "	V_next[comp] = V[comp] + ( dt * ( I_internal ) / C[comp] )%s;\n", Vnext_suffix.c_str()); ccde += tmps;
				// TODO explain the correction by showing how linearisation works
				const std::string CorrectI_suffix = Convert::Suffix(
					(Scales<Voltage>::native * Scales<Conductance>::native)
					.to(Scales<Current>::native)
				);
				sprintf(tmps, "		V_next[comp] = V[comp] + ( dt * ( I_internal + ( (V[comp] * G_internal)%s ) ) / C[comp] )%s;\n", CorrectI_suffix.c_str(), Vnext_suffix.c_str()); ccde += tmps;
				// TODO: NEURON runs all currents twice to approximate dI/dV;
				// perhaps do something like that and add the resulting diagonal to RC constant of compartments
			}
			
			ccde += "	}";
			
			// ccde += tab+"}\n";
			
			return true;
		};
		// then possibly integrate the cable equation as a separate step
		auto ImplementPostInternalCableEqIntegration = [ &config ](
			const SignatureAppender_Table &AppendMulti,
			const std::string &for_what,
			const std::string &tab,
			const SimulatorConfig::CableEquationSolver &cell_cable_solver,
			CellInternalSignature::PhysicalCell::CableSolverImplementation &cabl_impl,
			std::string &code
		){
			char tmps[2000];
			// now add a whole-cell solver here
			// TODO inline the process for small cells
			if( cell_cable_solver == SimulatorConfig::CABLE_BWD_EULER ){
				cabl_impl.Index_BwdEuler_OrderList   = AppendMulti.ConstI64("Bwd Euler Elimination Order");
				cabl_impl.Index_BwdEuler_ParentList = AppendMulti.ConstI64("Bwd Euler Elimination Parent");
				// work diagnonal and inv rc constant have already been ellocated and exposed
				
				// NB: if initial_state, then Vcomp_next := Vcomp through compartment's internal dynamics
				code += tab+"if(!initial_state){\n";
				// code += tab+"	printf(\"diagonal!! \\n\");\n";
				
				sprintf(tmps, "	const long long Compartments = cell_state_table_f32_sizes[%zd]; //same for all parallel arrays\n", cabl_impl.Index_BwdEuler_WorkDiagonal ); code += tmps;
				sprintf(tmps, "	const Table_I64 Order  = cell_const_table_i64_arrays[%zd];\n", cabl_impl.Index_BwdEuler_OrderList); code += tmps;
				sprintf(tmps, "	const Table_I64 Parent = cell_const_table_i64_arrays[%zd];\n", cabl_impl.Index_BwdEuler_ParentList); code += tmps;
				code += "	Table_F32 D = PerComp_InvRC_Diagonal;";
				// code += tab+"	printf(\"diagonal!!! \\n\");\n";	
				
				// diagonal = [ 1 - (1./ Capacitances[i]) * Y[i][i] * dt for i in range(Cells) ]
				// X = np.array(Voltages)
				const std::string Rate_suffix = Convert::Suffix(
					( Scales<Frequency>::native * Scales<Time>::native ) // to unitless
				);
				const std::string RCT_suffix = Convert::Suffix(
					( Scales<Time>::native / ( Scales<Resistance>::native * Scales<Capacitance>::native ) ) // to unitless
				);
				code += tab+"for(long long comp_seq = 0; comp_seq < Compartments; comp_seq++){\n";
				// sprintf(tmps, "		D[comp_seq] = 1 + DperT[comp_seq] * dt %s;\n", Rate_suffix.c_str() ); code += tmps;
				// sprintf(tmps, "		D[comp_seq] = 1 + PerComp_InvRC_Axial[comp_seq] * dt %s;\n", Rate_suffix.c_str() ); code += tmps;
				
				// convert the collected diagonal 1/RC to the unitless diagonal of the timestep matrix 
				sprintf(tmps, "		D[comp_seq] = 1 + D[comp_seq] * dt %s;\n", Rate_suffix.c_str() ); code += tmps;
				code += tab+"}\n";
				
				// lo_diagonal = [ -(1./ Capacitances[j]) * Y[i][j] * dt  for i in range(Cells) for j in (lastbranch[i],) ]
				// hi_diagonal = [ -(1./ Capacitances[i]) * Y[j][i] * dt  for i in range(Cells) for j in (lastbranch[i],) ]
				// code += tab+"	printf(\"diagonal#! \\n\");\n";	
				
				// Forward elimination
				// for ordi in range(Cells - 1):
				// 	i = elimination_order[ordi]
				// 	j = lastbranch[i]
				// 	# print "Eliminate %d %d" % ( i, j )
				// 	ratio = lo_diagonal[i] / diagonal[i] 
				// 	diagonal[j] -= ratio * hi_diagonal[i]
				// 	lo_diagonal[i] = 0
				// 	X[j] -= ratio * X[i]
				code += tab+"for( long long comp_seq = 0; comp_seq < Compartments - 1; comp_seq++ ){\n";
				code += tab+"	long long i = Order[comp_seq];\n";			
				code += tab+"	long long j = Parent[i];\n";
				code += tab+"	long long idx = ( ( i > j ) ? i : j );\n";
				code += tab+"	float R = R_Axial[idx];\n";
				sprintf(tmps, "		float Ui = - dt/( R * C[i]) %s;\n", RCT_suffix.c_str() ); code += tmps;
				sprintf(tmps, "		float Uj = - dt/( R * C[j]) %s;\n", RCT_suffix.c_str() ); code += tmps;
				code += tab+"	float Li = Uj;\n"; // really, check the math
				code += tab+"	float ratio = Li/D[i];\n";
				code += tab+"	D[j] -= ratio * Ui;\n";
				code += tab+"	V_next[j] -= ratio * V_next[i];\n";
				if( config.debug ){
				code += tab+"	printf(\"%lld %lld %g %g \\n\", i, j, D[i], V_next[i]);\n";	
				}		
				code += tab+"}\n";
				
				// print "After forward elimination:"
				// print diagonal
				// print lo_diagonal
				// print hi_diagonal
				// print X
				
				// i = elimination_order[-1]
				// X[i] = X[i] / diagonal[i]
				// # print "Last element", i, X[i]
				// diagonal[i] = 1
				code += tab+"long long i = Order[ Compartments - 1 ];\n";			
				code += tab+"V_next[i] = V_next[i] / D[i];\n";
				
				// Back-substitution
				// for ordi in range(Cells-2, -1, -1):
				// 	i = elimination_order[ordi]
				// 	j = lastbranch[i]
				// 	# print "Substitute %d %d" % ( i, j )
				// 	X[i] = (X[i] - hi_diagonal[i] * X[j] ) / diagonal[i]
				// 	hi_diagonal[i] = 0
				// 	diagonal[i] = 1
				code += tab+"for( long long comp_seq = Compartments - 2; comp_seq >= 0 ; comp_seq-- ){\n";
				code += tab+"	long long i = Order[comp_seq];\n";			
				code += tab+"	long long j = Parent[i];\n";
				code += tab+"	long long idx = ( ( i > j ) ? i : j );\n";
				code += tab+"	float R = R_Axial[idx];\n";
				sprintf(tmps, "		float Ui = - dt/( R * C[i]) %s;\n", RCT_suffix.c_str() ); code += tmps;
				code += tab+"	V_next[i] = ( V_next[i] - Ui * V_next[j] ) / D[i];\n";
				if( config.debug ){
				code += tab+"	printf(\"%lld %lld %g \\n\", i, j, V_next[i]);\n";	
				}
				code += tab+"}\n";
				
				code += tab+"}\n";
			}
			else{
				// no post-internal solver needed
			}
			
			return true;	
		};
		// lastly, add post-integration event handlers
		auto AllocateCreatePostIntegrationCode = [ &config, &ImplementSpikeSender, &cell_seq ](
			const SignatureAppender_Table &AppendMulti,
			size_t comp_seq,
			const std::string &for_what,
			const CellInternalSignature::PhysicalCell::CompartmentDefinition &comp_def,
			CellInternalSignature::PhysicalCell::CompartmentImplementation &comp_impl,
			std::string &code
		){
			(void) config; // just in case
			//create a vector for possible distribution of spike messages from presynaptic compartments
			if( comp_def.spike_output ){
				
				// check for Vth here too
				if( !std::isfinite(comp_def.Vt) ){
					// TODO more descriptive, with names, proper ids etc. mapping to segment especially
					printf("error: Cell type %zd compartment %zd has undefined Vthreshold, cannot use as spike source!\n", cell_seq, comp_seq);
					// TODO put check higher up, in NeuroML parsing (but will it hurt parsing efficiency? probably not?)
					return false;
				}
				
				if( !ImplementSpikeSender(
					"V[comp] <  V_threshold[comp] && V_threshold[comp] < V_next[comp]",
					AppendMulti,
					for_what,
					comp_impl.spiker, code
				) ) return false;
			}
			
			// done
			return true;
		};
		
		auto &compartment_grouping = pig.compartment_grouping = CellInternalSignature::CompartmentGrouping::AUTO;
		if( compartment_grouping == CellInternalSignature::CompartmentGrouping::AUTO ){
			// TODO more sophisticated analysis, cmd line/config options, etc etc.
			if( number_of_compartments <= 10 ) compartment_grouping = CellInternalSignature::CompartmentGrouping::FLAT;
			else compartment_grouping = CellInternalSignature::CompartmentGrouping::GROUPED;
		}
		
		// when compartments are fully flattened
		if( compartment_grouping == CellInternalSignature::CompartmentGrouping::FLAT ){
			sig.code += ExposeSubitemContext("local", "global", "\t");
		
			// add integration of dynamics
			for( int comp_seq = 0; comp_seq < number_of_compartments; comp_seq++ ){
				sprintf(tmps, "	// Internal Code for compartment %d\n", (int)comp_seq); sig.code += tmps;
				
				sprintf(tmps, "	{ int comp = %d;\n", (int)comp_seq); sig.code += tmps;
				
				std::string for_what = "Comp "+ itos(comp_seq);
				
				std::string intracomp_code;
				if( !ImplementInternalCompartmentIntegration(
					AppendSingle_CellScope, AppendMulti_CellScope, DescribeLemsInline_CellScope,
					for_what, tab,
					true,
					cell_cable_solver, bioph,
					comp_definitions[comp_seq], pig.comp_implementations[comp_seq], sig.cell_wig.random_call_counter,
					intracomp_code
				) ) return false;
				sig.code += intracomp_code;
				sig.code += "}";
				
				sprintf(tmps, "	// Internal Code for compartment %d end\n", (int)comp_seq); sig.code += tmps;
			}
			
			std::string cable_solver_code;
			if( !ImplementPostInternalCableEqIntegration(
				AppendMulti_CellScope, "", tab, cell_cable_solver,
				pig.cable_solver_implementation, cable_solver_code
			) ) return false;
			sig.code += cable_solver_code;
			
			for( int comp_seq = 0; comp_seq < number_of_compartments; comp_seq++ ){
				sprintf(tmps, "	// PostUpdate Code for compartment %d\n", (int)comp_seq); sig.code += tmps;
				sprintf(tmps, "	{ int comp = %d;\n", (int)comp_seq); sig.code += tmps;
				
				std::string for_what = "Comp "+ itos(comp_seq);
				
				std::string post_code;
				if( !AllocateCreatePostIntegrationCode( AppendMulti_CellScope, comp_seq, for_what, comp_definitions[comp_seq], pig.comp_implementations[comp_seq], post_code ) ) return false;
				sig.code += post_code;
				
				sprintf(tmps, "\t}\n\t// PostUpdate Code for compartment %d end\n", (int)comp_seq); sig.code += tmps;
			}
			
		}
		else if( compartment_grouping == CellInternalSignature::CompartmentGrouping::GROUPED ){
			
			auto &gp = pig.comp_group_impl;
			
			// compartment signature analysis
			// TODO set analysis and decision as a preliminary stage
			
			// well, a structural signature may be confused by having different mechanisms that are functionally equivalent (like different exp synapse types)
			// These may be generated by randomized instances of models, for example; so it's not so much up to the author to de-duplicate them
			// To confidently assert equivalence, all sub-components etc. of the mechanisms must be compared
			// So instead of peeforming the symbolic comparison, one can just
			// TODO find an elegant way to avoid false negative equivalence confusions LATER;
			// for example, if the LEMS components are the same in practice (just with different parms),  
			// it should work fine without too much effort
			
			// merging different compartment types with conditional enabling LATER
			// 	due to the pitfalls described above
			
			// isolate/generate the per-compartment code block, to group identical ones
			auto AllocateCreateFullCompartmentCode = [
				&ImplementInternalCompartmentIntegration, &AllocateCreatePostIntegrationCode,
				&model, &cell_cable_solver, &bioph
			](
				size_t comp_seq,
				const std::string &for_what,
				const std::string &tab,
				const CellInternalSignature::CompartmentDefinition &comp_def,
				CellInternalSignature::CompartmentImplementation &comp_impl,
				CellInternalSignature::WorkItemDataSignature &wig,
				std::string &intracomp_code, std::string &post_code
			){
				
				SignatureAppender_Single AppendSingle_CompScope( wig );
				SignatureAppender_Table AppendMulti_CompScope( wig );
				
				InlineLems_AllocatorCoder DescribeLemsInline_CompScope( model, wig.random_call_counter, AppendSingle_CompScope, AppendMulti_CompScope );
				
				if( !ImplementInternalCompartmentIntegration(
					AppendSingle_CompScope, AppendMulti_CompScope, DescribeLemsInline_CompScope,
					for_what, tab,
					false,
					cell_cable_solver, bioph,
					comp_def, comp_impl, wig.random_call_counter,
					intracomp_code
				) ) return false;
				
				if( !AllocateCreatePostIntegrationCode( AppendMulti_CompScope, comp_seq, for_what, comp_def, comp_impl, post_code ) ) return false;
				
				return true;
			};
			
			gp.distinct_compartment_types.clear();
			
			std::unordered_map< std::string, Int > compartment_code_hash_table;
			for( int comp_seq = 0; comp_seq < number_of_compartments; comp_seq++ ){
				// printf("comp %zd\n", comp_seq);
				// useful dummies
				CellInternalSignature::WorkItemDataSignature wig;
				CellInternalSignature::CompartmentImplementation comp_impl;
				// only this part matters
				std::string intracomp_code, post_code;
				
				if(!( AllocateCreateFullCompartmentCode(
					comp_seq,
					"", tab,
					comp_definitions[comp_seq],
					comp_impl,
					wig,
					intracomp_code, post_code
				) )) return false;
				
				auto key = intracomp_code + post_code;
				
				if( compartment_code_hash_table.count(key) > 0 ){
					gp.distinct_compartment_types.at( compartment_code_hash_table.at(key) ).Addd(comp_seq);
				}
				else{
					Int new_idx = (Int) gp.distinct_compartment_types.size();
					IdListRle l;
					l.Addd(comp_seq);
					gp.distinct_compartment_types.push_back(l);
					compartment_code_hash_table[key] = new_idx;
				}
				// printf("comp %zd end\n", comp_seq);
			}
			
			if(config.verbose){
			printf("Compartment types:\n");
			for( auto complist : gp.distinct_compartment_types ){
				printf( "\t%s\n", complist.Stringify().c_str() );
			}
			}

			// allocate the tables
			
			gp.Index_Coff    = AppendMulti_CellScope.ConstI64("Compartment Scalar CF32 Offset");
			gp.Index_Soff    = AppendMulti_CellScope.ConstI64("Compartment Scalar SF32 Offset");
			gp.Index_CF32off = AppendMulti_CellScope.ConstI64("Compartment Table  CF32 Offset");
			gp.Index_SF32off = AppendMulti_CellScope.ConstI64("Compartment Table  SF32 Offset");
			gp.Index_CI64off = AppendMulti_CellScope.ConstI64("Compartment Table  CI64 Offset");
			gp.Index_SI64off = AppendMulti_CellScope.ConstI64("Compartment Table  SI64 Offset");
			gp.Index_Roff    = AppendMulti_CellScope.ConstI64("Compartment RNG Offset");
			
			sig.code +=	tab+"const Table_I64 Comp_Coff    = cell_const_table_i64_arrays["+itos( gp.Index_Coff    )+"];\n";
			sig.code +=	tab+"const Table_I64 Comp_Soff    = cell_const_table_i64_arrays["+itos( gp.Index_Soff    )+"];\n";
			sig.code +=	tab+"const Table_I64 Comp_CF32off = cell_const_table_i64_arrays["+itos( gp.Index_CF32off )+"];\n";
			sig.code +=	tab+"const Table_I64 Comp_SF32off = cell_const_table_i64_arrays["+itos( gp.Index_SF32off )+"];\n";
			sig.code +=	tab+"const Table_I64 Comp_CI64off = cell_const_table_i64_arrays["+itos( gp.Index_CI64off )+"];\n";
			sig.code +=	tab+"const Table_I64 Comp_SI64off = cell_const_table_i64_arrays["+itos( gp.Index_SI64off )+"];\n";
			sig.code +=	tab+"const Table_I64 Comp_Roff    = cell_const_table_i64_arrays["+itos( gp.Index_Roff    )+"];\n";
			
			gp.preupdate_codes.resize ( gp.distinct_compartment_types.size() );
			gp.postupdate_codes.resize( gp.distinct_compartment_types.size() );
			gp.Index_CompList.resize  ( gp.distinct_compartment_types.size() );
			
			// TODO refactor better, restrict capture
			auto LoopOverCompartmentsCode = [ & ](
				size_t comptype_seq,
				std::string inner_code,
				std::string &ctde
			){
				// interlace it with compartment's tables, why not
				const auto Index_List = gp.Index_CompList[ comptype_seq ];
				
				ctde += tab+"// Internal Code for compartment type "+itos(comptype_seq)+"\n";
				ctde +=	tab+"{\n";
				
				ctde +=	tab+"const Table_I64 Comp_List    = cell_const_table_i64_arrays["+itos( Index_List    )+"];\n";
				ctde +=	tab+"const long long Type_Compartments    = cell_const_table_i64_sizes ["+itos( Index_List    )+"];\n";
				
				ctde +=	tab+"for( long long CompIdx = 0; CompIdx < Type_Compartments; CompIdx++ ){\n";
				
				ctde +=	tab+"	int comp = (int) Comp_List[CompIdx];\n";
				
				
				ctde += tab+"	const long long const_comp_index      = Comp_Coff   [comp];\n";
				ctde += tab+"	const long long state_comp_index      = Comp_Soff   [comp];\n";
				ctde += tab+"	const long long table_cf32_comp_index = Comp_CF32off[comp];\n";
				ctde += tab+"	const long long table_ci64_comp_index = Comp_CI64off[comp];\n";
				ctde += tab+"	const long long table_sf32_comp_index = Comp_SF32off[comp];\n";
				ctde += tab+"	const long long table_si64_comp_index = Comp_SI64off[comp];\n";
				
				ctde += tab+"	const long long rng_offset            = Comp_Roff   [comp];\n";
				ctde += tab+"	\n";
				
				ctde += ExposeSubitemContext("comp", "cell", "\t");
				ctde += CloneSubitemIndices("local", "comp", "\t");
				
				ctde += ExposeSubitemContext("local", "cell", "\t");
				
				ctde += inner_code;
				
				ctde +=	tab+"}\n";
				
				
				ctde +=	tab+"}\n";
				ctde+= tab+"// Internal Code for compartment type "+itos(comptype_seq)+" end\n";
				
				return true;
			};
			
			
			// keep tne necessary offsets for deduplication (localization) of each compartment's internals
			// NB: nComps == number_of_compartments here, but it could break off in case only a part of the cell is being handled this way...
			int nComps = (int)pig.comp_implementations.size();
			gp.r_off   .resize( nComps );
			gp.c_off   .resize( nComps );
			gp.s_off   .resize( nComps );
			gp.cf32_off.resize( nComps );
			gp.sf32_off.resize( nComps );
			gp.ci64_off.resize( nComps );
			gp.si64_off.resize( nComps );
			
			for( size_t comptype_seq = 0; comptype_seq < gp.distinct_compartment_types.size(); comptype_seq++ ){
				
				gp.Index_CompList[ comptype_seq ] = AppendMulti_CellScope.ConstI64("List of Type "+itos(comptype_seq)+" Compartments");
				
				bool first_compartment = true;
				for( int comp_seq : gp.distinct_compartment_types[comptype_seq].toArray() ){
					// first allocate the wigs and full codes for every compartment
					// superfluous codes may be dropped for similar compartments
					
					std::string for_what = "Comp "+ itos(comp_seq);
					
					std::string intracomp_code, post_code;
					
					// NB: when creating compartment implementations, the implementation is defined based on the wig it's appended to.
					// So, on one side fictitious wigs are useful to deduplicate code depending on the implementation's indices,
					// on the oher side the actual implementations need to be stored in an absolute manner, so pigs can be accessed in the same way whether flat or grouped.
					// 	(The other option is to do the necessary translation in each and every pig-accessing interface)
					// BEWARE when implementing vectors that refer to work-item indices,
					// for the illusion must be kept inside the deduplicated code, and removed outside this part of the kernel.
					// TODO add a more elegant formulation.
					
					// Use a fictitious, standalone (ie local namespace) wig, just for the code
					if( first_compartment ){
						
						// it is like a wig, but fake-er
						CellInternalSignature::WorkItemDataSignature faux_wig;
						// could also overwrite, beware of operations lacking overwrite semantics like map::insert !
						CellInternalSignature::PhysicalCell::CompartmentImplementation faux_comp_impl;
						
						std::string for_what = "Comp Type "+ itos(comptype_seq); // this is the code that will be used within the per compartment type loop
						
						if(!( AllocateCreateFullCompartmentCode(
							comp_seq,
							for_what, tab,
							comp_definitions[comp_seq],
							faux_comp_impl, // pig.comp_implementations[comp_seq],
							faux_wig,
							intracomp_code, post_code
						) )) return false;
						
						
						gp.preupdate_codes[comptype_seq] = intracomp_code;
						gp.postupdate_codes[comptype_seq] = post_code;
					}
					
					// but actually allocate the implementation, by appending on the cell sig
					
					auto &cell_wig = sig.cell_wig;
					
					gp.r_off   [comp_seq] = ( cell_wig.random_call_counter   );
					gp.c_off   [comp_seq] = ( cell_wig.constants       .size() );
					gp.s_off   [comp_seq] = ( cell_wig.state           .size() );
					gp.cf32_off[comp_seq] = ( cell_wig.tables_const_f32.size() );
					gp.sf32_off[comp_seq] = ( cell_wig.tables_state_f32.size() );
					gp.ci64_off[comp_seq] = ( cell_wig.tables_const_i64.size() );
					gp.si64_off[comp_seq] = ( cell_wig.tables_state_i64.size() );
					
					// could decouple code generation from data signature allocation, for efficiency LATER
					// In theory, a clever enough compiler could flatten the unnecessary code generation out of existence,
					// since it has no effect (other than allocating memory)
					if(!( AllocateCreateFullCompartmentCode(
						comp_seq,
						for_what, tab,
						comp_definitions[comp_seq],
						pig.comp_implementations[comp_seq],
						cell_wig,
						intracomp_code, post_code
					) )) return false;
					
					first_compartment = false;
				}
			}
			
			// and append the internal dynamics codes
			for( size_t comptype_seq = 0; comptype_seq < gp.distinct_compartment_types.size(); comptype_seq++ ){
				
				std::string comptype_inner_code;
				if( !LoopOverCompartmentsCode( comptype_seq, gp.preupdate_codes[comptype_seq], comptype_inner_code ) ) return false;
				
				sig.code += comptype_inner_code;
			}
			// and the cable solver, if it's a separate step
			std::string cable_solver_code;
			if( !ImplementPostInternalCableEqIntegration(
				AppendMulti_CellScope, "", tab, cell_cable_solver,
				pig.cable_solver_implementation, cable_solver_code
			) ) return false;
			sig.code += cable_solver_code;
			
			// and finally the postupdate codes
			for( size_t comptype_seq = 0; comptype_seq < gp.distinct_compartment_types.size(); comptype_seq++ ){
				
				sig.code += tab+"// PostUpdate Code for compartment type "+itos(comptype_seq)+"\n";
				
				std::string comptype_outer_code;
				if( !LoopOverCompartmentsCode( comptype_seq, gp.postupdate_codes[comptype_seq], comptype_outer_code ) ) return false;
				
				sig.code += comptype_outer_code;
			}
			
			// done for now, copy vectors when initializing later on
		}
		else{
			printf("internal error: unknown compartment grouping %d for cell type %d", (int) compartment_grouping, (int) cell_seq);
			return false;
		}
		EmitWorkItemRoutineFooter( sig.code );
		EmitKernelFileFooter( sig.code );
		// done with code generation for this work item
		
		}
		else if( cell_type.type == CellType::ARTIFICIAL ){
			char tmps[10000]; // buffer for a single code line
			
			const auto &cell = cell_type.artificial;
			
			auto &cell_wig = sig.cell_wig;
			auto &aig = sig.artificial_cell;
			
			auto &ccde = sig.code;
			
			
			auto &AppendSingle = AppendSingle_CellScope;
			auto &AppendMulti  = AppendMulti_CellScope;
			auto &DescribeLemsInline  = DescribeLemsInline_CellScope;
			
			printf("Generating code for %s...:\n", sig.name.c_str());
			
			EmitKernelFileHeader( sig.code );
			EmitWorkItemRoutineHeader( sig.code );
			
			const std::string tab = "\t";			
			std::string for_what = "Cell";		
			
			// expose the local context of the artificial cell
			sig.code += ExposeSubitemContext("local", "global", "\t");
			
			// per cell RNG
			ImplementRngSeed(
				AppendSingle_CellScope,
				"", tab,
				"local",
				sig.common_in_cell.cell_rng_seed,
				sig.code
			);	
			sig.code += "	const int rng_object_id = cell_rng_seed;\n";
			
			ccde   += tab+"char spike_in_flag = 0;\n";
			ccde   += tab+"char spike_out_flag = 0;\n";
			
			auto MaybeDetermineComponentVoltage_Lems = [ &model ](
				const ComponentInstance &comp_inst,
				CellInternalSignature::ArtificialCell &aig
			){
				const auto &comp_type = model.component_types.get(comp_inst.id_seq);
				aig.Index_Statevar_Voltage = -1;
				const auto &voltage_thing_seq = comp_type.common_exposures.membrane_voltage;
				// printf("hello %zd %zd\n", comp_inst.id_seq, voltage_thing_seq);
				if( voltage_thing_seq >= 0 ){
					const auto &voltage_thing = comp_type.exposures.get(voltage_thing_seq);
					// printf("heo %zd %d, %s %zd\n", voltage_thing.seq, voltage_thing.type, voltage_thing.Stringify().c_str(), aig.component.statevars_to_states.size());
					
					if( voltage_thing.type == ComponentType::Exposure::STATE ){
						aig.Index_Statevar_Voltage = aig.component.statevars_to_states[ voltage_thing.seq ].index;
					}
					else{
						// printf("artificial cell noes not expose voltage as state variable\n");
						// return false; // perhaps reconsider LATER
					}
				}
			};
			auto MaybeDetermineComponentVoltage = [ &MaybeDetermineComponentVoltage_Lems ](
				const ArtificialCell &cell,
				CellInternalSignature::ArtificialCell &aig
			){
				if( !cell.component.ok() ){
					// expose voltage of handmade cells LATER
					return;
				}
				else{
					// LEMS implementation
					
					MaybeDetermineComponentVoltage_Lems(cell.component, aig);
					return;
				}
			};
			
			
			if( cell.type == ArtificialCell::SPIKE_SOURCE ){
				
				// Now things are about to get ugly, because inputs are typically implemented as tables (and being able to specialize is important for performance)
				// At the same time, having inlined artificial cells is also important for performance.
				// So this either takes developing one set of completely layout-agnostic implementations what still do not sacrifice performance,
				// or two sets of implemetations, keeping one for the "spike source as artificial cell" case.
				// The second path is chosen, since spike sources don't really make sense as input sources,
				//   and they make much more sense as standalone artificial cells
				
				// internal dynamics only
				// Vcomp is not needed since it's a lonely input all by itself, not attached to any membrane 
				
				const InputSource &input = input_sources.get(cell.spike_source_seq);
				if( input.type == InputSource::SPIKE_LIST ){
					const auto &for_that = for_what;
					std::string for_what = for_that + " Spike List";
					
					auto &inpimpl = aig.inpimpl;
					char tmps[1000];
					
					// for each instance of the same input source type:
					// a slice of a common spike time vector, and a start index to begin from for each instance
					// and a state variable to which index is coming up for each instance
					size_t table_Times = inpimpl.Table_SpikeListTimes  = AppendMulti.Constant(for_what+" Spike Times");
					size_t table_Posit = inpimpl.Table_SpikeListPos    = AppendSingle.StateVariable( 0, for_what+" Spike Index Position Integer"); // should be an integer state, oh well
					
					// positions are initialized at initial tables time, yay!
		
					sprintf(tmps, "	const long long Instances = local_state_table_i64_sizes[%zd]; //same for all parallel arrays\n", inpimpl.Table_SpikeListPos ); ccde += tab+tmps;
					
					sprintf(tmps, "const float *Spike_Times = local_const_table_f32_arrays[%zd];\n", table_Times); ccde += tab+tmps;
					sprintf(tmps, "const float *Position  = &local_state    [%zd];\n", table_Posit); ccde += tab+tmps;
					sprintf(tmps, "      float *PositNext = &local_stateNext[%zd];\n", table_Posit); ccde += tab+tmps;
					
					// TODO wrap into a reqstring ?
					ccde   += tab+"{\n";
					
					bool safe_cast = true;
					// TODO profile, and do something both fast and standard
					if( safe_cast ){
						ccde   += tab+"int pos = (int) *Position;\n"; // safe, probably slow
					}
					else{
						// funny business may (in theory) lead to a trap value! beware!
						ccde   += tab+"union TypePun{ int i32; float f32; } cast;\n";
						ccde   += tab+"{ char static_assert[ sizeof(int) == sizeof(float) ]; }\n";
						ccde   += tab+"cast.f32 = *Position; int pos = cast.i32;\n";
					}
					
					ccde   += tab+"if( !initial_state ){\n";
						
					ccde   += tab+"	while( time_f32 >= Spike_Times[pos] ){\n";
					ccde   += tab+"		spike_out_flag |= 1;\n";
					ccde   += tab+"		pos++;\n";
					ccde   += tab+"	}\n";
						
					ccde   += tab+"}\n";
					ccde   += tab+"else{\n";
					ccde   += tab+"	pos = 0; // initialize\n";
					ccde   += tab+"}\n";
					
					if( safe_cast ){
						ccde   += tab+"*PositNext = (float)pos;\n"; // safe, probably slow
					}
					else{
						ccde   += tab+"cast.i32 = pos; *PositNext = cast.f32;\n"; // fast, probably unsafe 
					}
					
					ccde   += tab+"}\n";
					
				}
				else{
					// no other native, only LEMSified implementations for now
					if( input.component.ok() ){
						const ComponentInstance &comp_inst = input.component;
						ccde += tab+"{\n";
				
						ccde += DescribeLemsInline.SingleInstance( comp_inst, tab, for_what, aig.component, config.debug );
						
						ccde += tab+"spike_out_flag |= Lems_eventout_spike;\n";
						
						ccde += tab+"}\n";
					}
					else{
						printf("Unknown native (input as artificial cell) type\n");
						return false;
					}
				}
				
			}
			else{
				// must allocate first, for voltage to be exposed
				if( cell.component.ok() ){
					const auto &compinst = cell.component;
					const auto &comptype = model.component_types.get(compinst.id_seq);
					
					aig.component = DescribeLems::AllocateSignature(comptype, compinst, &AppendSingle, for_what + " LEMS");
					
				}
				else{
					printf("Unknown native artificial cell type\n");
					return false;
				}
				// must expose voltage;
				// in fact, syn current should be computed:
				// 	- after the LEMS vars it depends on 
				// 	- , and *before* the LEMS vars that depend on it
				// 	good luck with that
				// for now, just expose voltage which is sure (?) to be available
				// TODO something more elegant
				MaybeDetermineComponentVoltage(cell, aig);
				
				if( aig.Index_Statevar_Voltage >= 0 ){
					sprintf(tmps, "	const float Vcomp = local_state[%zd]; \n", (size_t)aig.Index_Statevar_Voltage); ccde += tmps;
				}
				
				// add synapses
				ccde += "	// Current from synapses\n";
				ccde += "	float I_synapses_total = 0;\n";
				ccde += "	float G_synapses_total = 0;\n";
				for( const auto &keyval : synaptic_component_types_per_cell_per_compartment[cell_seq] ){
					for( const auto &id_id : keyval.second.toArray() ){
						if( !ImplementSynapseType( AppendSingle, AppendMulti, DescribeLemsInline, cell_wig.random_call_counter, for_what + " Synapse type "+std::to_string(id_id), tab, id_id, aig.synapse, ccde ) ) return false;
					}
				}
				
				// add inputs
				ccde += "	// Current from inputs\n";
				ccde += "	float I_input_total = 0;\n";
				ccde += "	float G_input_total = 0;\n";
					
				// generate tables and code for each input type
				for( const auto &keyval : input_types_per_cell_per_compartment[cell_seq] ){
					for( const auto &id_id : keyval.second.toArray() ){			
						if( !ImplementInputSource( AppendSingle, AppendMulti, DescribeLemsInline, cell_wig.random_call_counter, for_what + " Input type "+std::to_string(id_id), tab, id_id, aig.input, ccde ) ) return false;
					}
				}
				
				ccde += "	float external_current = I_synapses_total + I_input_total;\n";
				// unfortunately G_input_total can't be used to improve the integrator since arbitrary I dynamics are at play. (Unless the LEMS component exposes the interation of I into V rule to be accessible somehow...)
				
				// internal dynamics
				
				// no other native, only LEMSified implementations for now
				if( cell.component.ok() ){
					const auto &compinst = cell.component;
					const auto &comptype = model.component_types.get(compinst.id_seq);
					
					auto &component = aig.component;
					
					ccde += tab+"// LEMS assigned\n";
					std::string lemscode = DescribeLems::Assigned(comptype, model.dimensions, component, &AppendSingle_CellScope, for_what, tab, cell_wig.random_call_counter, config.debug && 0 );
					ccde += lemscode;
					
					// also add integration code here, to finish with component code (and get event outputs !)
					ccde += tab+"// integrate inline\n";
					std::string lemsupdate = DescribeLems::Update(comptype, model.dimensions, component, &AppendSingle, for_what, tab, cell_wig.random_call_counter, config.debug && 0 );
					ccde += lemsupdate;
					
					ccde += tab+"// expose inline\n";
					ccde += DescribeLems::Exposures( comptype, for_what, tab, config.debug && 0 );
					
					// also pass the spike output through, if there is one
					if( comptype.common_event_outputs.spike_out >= 0 ){
						ccde += tab+"spike_out_flag |= Lems_eventout_spike;\n";
					}
					
					// if( !ImplementInternalLemsDynamics( cell.component, aig, ccde ) ) return false;
				}
				else{
					printf("Unknown native artificial cell type\n");
					return false;
				}
				
				// NB Index_Statevar_Voltage should have been set by this point, if available
			
			}
			// send (exposure)
			if( spiking_outputs_per_cell_per_compartment.at( cell_seq ).size() > 0 ){
				// could be added regardless, but add it only when needed for now
				if( !ImplementSpikeSender(
					"!!spike_out_flag",
					AppendMulti,
					for_what,
					aig.spiker, ccde
				) ) return false;
			}
			
			EmitWorkItemRoutineFooter( sig.code );
			EmitKernelFileFooter( sig.code );
		}
		
		
		// printf("%s", sig.code.c_str());
		// printf("\n");
		
		//and printout the whole signature
		auto PrintWorkItemSignature = []( const CellInternalSignature::WorkItemDataSignature &wig ){
			printf("Constants:\n");
			for(size_t i = 0; i < wig.constants.size(); i++){
				printf("\t%20g\t", wig.constants[i]);
				if(wig.constants_names.count(i)){
					printf("%s", wig.constants_names.at(i).c_str());
				}
				printf("\n");
			}
			
			printf("States:\n");
			for(size_t i = 0; i < wig.state.size(); i++){
				printf("%zd\t%20g\t",i, wig.state[i]);
				if(wig.state_names.count(i)){
					printf("%s", wig.state_names.at(i).c_str());
				}
				printf("\n");
			}

			printf("Tables:\n");
			auto PrintTabSig = []( const char *tabtype, const std::vector<CellInternalSignature::TableInfo> &tabsig ){
				for(size_t i = 0; i < tabsig.size() ; i++){
					const auto &inf = tabsig.at(i);
					printf("\t%s %3zd:\t %s\n", tabtype, i, inf.Description().c_str());
				}
			};
			PrintTabSig("CF32", wig.tables_const_f32);
			PrintTabSig("CI64", wig.tables_const_i64);
			PrintTabSig("SF32", wig.tables_state_f32);
			PrintTabSig("SI64", wig.tables_state_i64);
			
			printf("\n");
			printf("\n");
		};
		if(config.verbose){
		PrintWorkItemSignature(sig.cell_wig);
		}
		
		
		// output model code for all present cells TODO
		std::string code_id = sig.name + "_code";
		std::string code_filename = code_id+ ".gen.c";
		std::string dll_filename = code_id+ ".gen.so";
		FILE *fout = fopen(code_filename.c_str(), "w");
		if(!fout){
			perror(code_filename.c_str());
			return false;
		}
		if( !fprintf(fout, "%s", sig.code.c_str() ) ){
			perror(code_filename.c_str());
			return false;
		}
		fclose(fout);
		
		// build the code
		timeval compile_start, compile_end;
		gettimeofday(&compile_start, NULL);
		
		
		std::string basic_flags =
			" -std=c11 -Wall"
			" -Wno-attributes"
			" -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function";
		std::string dll_flags = " -shared -fpic"; //" -shared -fpic -nodefaultlibs";
		std::string optimization_flags = " -Ofast";
		// Optimization flags are quirky between CPU archs.
		// When using GCC, refer to: https://gcc.gnu.org/onlinedocs/gcc/Submodel-Options.html#Submodel-Options
		
		// for arch detection, refer to:
		// https://stackoverflow.com/questions/66249936/detecting-cpu-architecture-compile-time
		// and https://sourceforge.net/p/predef/wiki/Architectures/
		// also detecting arm on visual c++: https://stackoverflow.com/questions/37251625/detect-arm-64-in-preprocessor
		// TODO move arch aliases to Common.h
		
		#if defined(__i386__) || defined(__i386) || defined(_M_IX86) || \
			defined(__x86_64__) || defined(_M_X64)
		// the old reliable setting for the typical PC user base, mcpu is deprecated alias for mtune here
		optimization_flags += " -march=native -mtune=native";
		#elif defined(__mips)
		// also for mips for some reason? similarly old? -mcpu is not supported, so
		optimization_flags += " -march=native -mtune=native";
		#elif defined(__s390__) || defined(__s390x__)
		// also for s390. funny enough
		optimization_flags += " -march=native -mtune=native";
		#elif defined(__arm__) || defined(_M_ARM) || \
			defined(__aarch64__) || defined(_M_ARM64)
		
		// Now it gets interesting because of Apple Clang...
		// for GCC, refer to https://community.arm.com/arm-community-blogs/b/tools-software-ides-blog/posts/compiler-flags-across-architectures-march-mtune-and-mcpu
		// for Clang, refer to: can't see any answer here https://clang.llvm.org/docs/UsersManual.html
		// also: https://maskray.me/blog/2022-08-28-march-mcpu-mtune
		
		// reportedly, -mcpu=native is supported on clang since 2018, so it should be available on all setups
		// let's hope this works
		optimization_flags += " -mcpu=native";
		#else
		// some other arch, this has a chance of working, but risks failing
		// the true solution for compiler detection is some sort of autotools...
		optimization_flags += " -mcpu=native";
		// probably works for gcc on powerpc: https://gcc.gnu.org/onlinedocs/gcc/RS_002f6000-and-PowerPC-Options.html
		
		// will still fail on old versions/new platforms, oh well
		#endif
		std::string fastbuild_flags = " -O0";
		std::string asm_flags = " -S -masm=intel -fverbose-asm";
		std::string lm_flags = " -lm";
		if(config.use_icc){
			lm_flags = " -limf"; // TODO check if SVML should be added here too
		}
		if( !config.use_icc && config.tweak_lmvec ){
			// some GCC(glibc?) versions need this for vectorization https://sourceware.org/bugzilla/show_bug.cgi?id=20539
			// also, mvec must be linked first: https://sourceware.org/glibc/wiki/libmvec 
			lm_flags = " -lmvec -lm" ; // XXX must revert automatically, if compiler is old enough 
		}
		// TODO extra_flags, vec_report etc.
		
		std::string compiler_name = "icc";
		if(!config.use_icc){
			compiler_name = "gcc";
		}
		
		std::string code_quality_flags = optimization_flags;
		
		// don't bother with optimization if code is massive
		if( sig.code.size() > 1024 * 1024LL ){
			if(config.verbose){
				printf("Choosing fast build due to code size..\n");
			}
			code_quality_flags = fastbuild_flags;
		}
		
		// Check if compiler is present
		// TODO XXX hoist this check higher up, to avoid overhead !
		// TODO more branching to pick the method to check presence LATER, for more compilers
		if( system((compiler_name + " --version").c_str()) != 0 ){
			std::string complaint_line = "Could not invoke '"+compiler_name+"' compiler! Make sure it is installed, and available on PATH.";
			
			std::string more_commentary;
			// maybe mention that gcc is the default (or auto selected) LATER
			
			// Give some instructions to the astonished user, though the most complete instructions should really be in the manual (when that is written)
			if(config.use_icc){
				more_commentary = "Check the instructions on how to set up ICC at Intel's website:\n"
				"https://software.intel.com/content/www/us/en/develop/articles/intel-system-studio-download-and-install-intel-c-compiler.html"
				"\nand on setting PATH:\n"
				"https://software.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/compiler-setup/using-the-command-line/specifying-the-location-of-compiler-components.html";
			}
			else{
				// gcc by default
				#if defined _WIN32
				more_commentary = "If a compiler is not already installed, a build for GCC on Windows can be downloaded from:\n";
				
				#if INTPTR_MAX == INT32_MAX
				more_commentary += "https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/8.1.0/threads-posix/sjlj/i686-8.1.0-release-posix-sjlj-rt_v6-rev0.7z";
				#elif INTPTR_MAX == INT64_MAX
				more_commentary += "https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Wing4/Personal%20Builds/mingw-builds/8.1.0/threads-posix/seh/x86_64-8.1.0-release-posix-seh-rt_v6-rev0.7z";
				#endif
				// but still may have to detect non-PC architecture ... LATER

				more_commentary += "\nUnpack the file anywhere, and add the unpacked <path ...>\\bin directory to EDEN's PATH.";
				#elif defined __linux__
				more_commentary = "GCC is usually already installed on Linux setups. If it is not installed, refer to your distribution's documentation on how to install the essentials for building from source.";
				#elif defined(__APPLE__)
				more_commentary = "A GCC-compatible compiler cn be installed with the Command Line Developer Tools for Mac. Run the following command on the Terminal to install:\n";
				more_commentary += "xcode-select --install\n\n";
				more_commentary += "Alternatively, the compiler used by default, GCC, can be installed through Homebrew for Mac OS X:\n";
				more_commentary += "brew install gcc";
				more_commentary += "\nRefer to http://brew.sh on how to set up Homebrew. (It may already be installed, in order to install Python 3.)";
				#else
					
				#endif
			}
			// also note ways to set path
			#if defined _WIN32
			more_commentary += "If using the command line, PATH can be set as follows:\n"
			"path <path to compiler executable>;%PATH%\n"
			"eden.exe ..."; // maybe argv[0], whatever
			#elif defined __linux__
			more_commentary += "If using the command line, PATH can be set as follows:\n"
			"PATH=<path to compiler executable>:$PATH eden ...";
			#endif
			
			more_commentary += "If using Python, PATH can be set as follows:\n"
			"os.environ[\"PATH\"] = <path to compiler executable> + os.pathsep + os.environ[\"PATH\"]\n"
			"runEden(...)";
			
			fprintf(stderr, "%s\n", complaint_line.c_str());
			if( !more_commentary.empty() ){
				fprintf(stderr, "%s\n", more_commentary.c_str());
			}
			
			return false;
		}
		
		// NOTE -lm must be put last, after other obj files (like source code) have stated their dependencies on libm
		// further reading: https://eli.thegreenplace.net/2013/07/09/library-order-in-static-linking
		std::string cmdline =     compiler_name + " " + basic_flags + dll_flags + code_quality_flags + " -o " + dll_filename + " " + code_filename + lm_flags;
		if(config.verbose) printf("%s\n", cmdline.c_str());
		std::string cmdline_asm = compiler_name + " " + basic_flags + dll_flags + code_quality_flags + asm_flags + " " + code_filename + lm_flags;
		if(config.output_assembly){
			if( system(cmdline_asm.c_str()) != 0 ){
				fprintf(stderr, "Could not build %s assembly\n", dll_filename.c_str());
				return false;
			}
		}
		if( system(cmdline.c_str()) != 0 ){
			fprintf(stderr, "Could not build %s\n", dll_filename.c_str());
			return false;
		}
		
		// load the code
		std::string function_name = "doit";
		IterationCallback callback = NULL;
		
		#if defined (__linux__) || defined(__APPLE__)
		void *dll_handle = dlopen(("./"+dll_filename).c_str(), RTLD_NOW);
		if(!dll_handle){
			fprintf(stderr, "Error loading %s: %s\n", dll_filename.c_str(), dlerror());
			return false;
		}
		*(void**)(& callback ) = dlsym(dll_handle, function_name.c_str()); // C-style voodoo to make a "valid" cast
		if(!callback){
			fprintf(stderr, "Error loading %s symbol %s: %s\n", dll_filename.c_str(), function_name.c_str(), dlerror());
			dlclose(dll_handle);
			return false;
		}
		#endif
		#ifdef _WIN32
		// TODO normalize paths to place dll's somewhere else than cwd !
		// TODO Unicode support, with MultiByteToWideChar
		HMODULE dll_handle = LoadLibraryA((".\\"+dll_filename).c_str());
		if(!dll_handle){
			DWORD errCode = GetLastError();
			fprintf(stderr, "Error loading %s: %s\n", dll_filename.c_str(), DescribeErrorCode_Windows(errCode).c_str());
			return false;
		}
		*(void**)(& callback ) = (void*)GetProcAddress(dll_handle, function_name.c_str());
		if(!callback){
			DWORD errCode = GetLastError();
			fprintf(stderr, "Error loading %s symbol %s: %s\n", dll_filename.c_str(), function_name.c_str(), DescribeErrorCode_Windows(errCode).c_str());
			FreeLibrary(dll_handle);
			return false;
		}
		#endif
		
		if(!callback){
			// which is already guarded against in platform specific code.
			// the only reason this should happen is if the platform is not supported
			fprintf(stderr, "Error loading %s: %s\n", dll_filename.c_str(), "internal error");
			return false;
		}
		sig.callback = callback;
		// LATER keep a set of dynamic libraries loaded, to cleanup
		// though it's pointless in this sort of application
		gettimeofday(&compile_end, NULL);
		if(config.verbose) printf("Compiled and loaded %s in %.2lf seconds\n", code_id.c_str(), TimevalDeltaSec(compile_start, compile_end));
		
		
		// done with cell type's signature
	}
	// LATER further specialize cell types with synapse and input components INSIDE the per-cell code block, for better legibility, but how?
	
	// now realize the model
	
	//------------------> Generate the data structures
	
	// Instantiate cells, with internal working sets onto global tables and extension tables onto inner table space
	// Simplest layout : Array of Structures, lay out states and constants for every single cell, side by side
	// TODO structure-of-arrays transformation, at least over block-wise intervals
	// TODO instantiate class constants
	// TODO instantiate more global constants
	
	typedef ptrdiff_t work_t;
	
	
	// symbolic reference to some point, on some neuron, under NeuroML
	// perhaps move this to be more general, in header LATER
	struct PointOnCellLocator{
		Int population;
		Int cell_instance;
		Int segment;
		Real fractionAlong;
		
		// LATER avoid interpreting and expanding all populations on all nodes, perhaps offload it to a pre-processing tool ( like ncx, opencortex and such )
		bool operator<( const PointOnCellLocator &rhs ) const {
			if( population    < rhs.population    ) return true;
			if( population    > rhs.population    ) return false;
			if( cell_instance < rhs.cell_instance ) return true;
			if( cell_instance > rhs.cell_instance ) return false;
			if( segment       < rhs.segment       ) return true;
			if( segment       > rhs.segment       ) return false;
			if( fractionAlong < rhs.fractionAlong ) return true;
			if( fractionAlong > rhs.fractionAlong ) return false;
			
			return false;
		}
		
		std::string toPresentableString() const {
			std::string ret;
			ret +=
				"(pop " + presentable_string(population)
				+ ", cell " + presentable_string(cell_instance)
				+ ", seg " + presentable_string(segment)
				+ ", frac " + presentable_string(fractionAlong)
				+ ")";
			return ret;
		}
		
		void toEncodedString( std::string &out_str ) const {
			out_str +=
				accurate_string(population) + " "
				+ accurate_string(cell_instance) + " "
				+ accurate_string(segment) + " "
				+ accurate_string(fractionAlong)
			;
		}
		bool fromEncodedString( const char *in_str ){
			if( sscanf( in_str, "%ld %ld %ld %f", &population, &cell_instance, &segment, &fractionAlong ) != 4 ) return false;
			// would do more checking, if it wasn't internally used
			return true;
		}
	};
	
	#ifdef USE_MPI
	
	// get lists of what to be sent, in model-specific symbolic references -
	// not references to realized work items, because each node may handle this differently
	
	// symbolic reference to a DataWriter
	struct DawRef{
		// just its position on the NeuroML-supplied list
		Int daw_seq;
		Int col_seq;
		
		bool operator<( const DawRef &rhs ) const {
			if( daw_seq    < rhs.daw_seq    ) return true;
			if( daw_seq    > rhs.daw_seq    ) return false;
			if( col_seq    < rhs.col_seq    ) return true;
			if( col_seq    > rhs.col_seq    ) return false;
			
			return false;
		}
		
		std::string toPresentableString() const {
			std::string ret;
			ret +=
				"(daw " + presentable_string(daw_seq)
				+ ", col " + presentable_string(col_seq)
				+ ")";
			return ret;
		}
		
		void toEncodedString( std::string &out_str ) const {
			out_str +=
				accurate_string(daw_seq) + " "
				+ accurate_string(col_seq)
			;
		}
		bool fromEncodedString( const char *in_str ){
			if( sscanf( in_str, "%ld %ld ", &daw_seq, &col_seq ) != 2 ) return false;
			// would do more checking, if it wasn't internally used
			return true;
		}
		
	};
	
	
	// list of what this node sends to peers that need it
	struct SendList{
		// Order in vectors is same as order in sent packet
		// (and is the order that the receiver requested)
		std::vector<PointOnCellLocator> vpeer_sources;
		std::vector<DawRef> daw_refs;
		std::vector<PointOnCellLocator> spike_sources;
	};
	
	// list of what this node needs from other peers
	struct RecvList{
		
		// existing references to state variables, to be remapped to the mirror buffers according to type and point on cell
		std::map< PointOnCellLocator, std::vector< TabEntryRef_Packed > > vpeer_refs;
		
		// positions of trigger buffers, to be updated by spikes originating from points on cell
		std::map< PointOnCellLocator, std::vector< TabEntryRef_Packed > > spike_refs;
		
		// if a logging node, it needs to record values on remote work items
		std::set< DawRef > daw_refs;
	};
	
	// not the send/recv buffers themselves, but lists of references
	std::map< int, SendList > send_lists; // per adjacent node
	std::map< int, RecvList > recv_lists; // per adjacent node
	
	// put a provisional table entry, and keep track of these dependencies
	// such table entries _will_ be remapped to point to the corresponding positions on send/recv buffers
	auto AppendRemoteDependency_Vpeer = [ &tabs, &recv_lists ]( const PointOnCellLocator &loc, int remote_node, size_t glob_tab_Vpeer ){
		
		auto &table = tabs.global_tables_const_i64_arrays.at(glob_tab_Vpeer);
		
		long long entry = table.size();
		
		TabEntryRef_Packed packed_id = GetEncodedTableEntryId( glob_tab_Vpeer, entry );
		
		recv_lists[remote_node].vpeer_refs[loc].push_back(packed_id);	
		
		long long temp_id = -100 - remote_node;
		table.push_back(temp_id);
		
		return true;
	};
	auto AppendRemoteDependency_Spike = [ &recv_lists ]( const PointOnCellLocator &loc, int remote_node, TabEntryRef_Packed trig_buf_ref ){
		// printf("remotespike %d, %s, %llx, %zd\n", remote_node, loc.toPresentableString().c_str(), trig_buf_ref, recv_lists[remote_node].spike_refs[loc].size() );
		recv_lists[remote_node].spike_refs[loc].push_back(trig_buf_ref);	
		
		return true;
	};
	auto AppendRemoteDependency_DataWriter = [ &recv_lists ]( const DawRef &dawref, int remote_node ){
		
		recv_lists[remote_node].daw_refs.insert(dawref);	
		
		return true;
	};
	
	// like  workunit_per_cell_per_population, but explicitly referring to local node's fragment of the model
	std::vector< std::map<Int, work_t> >  local_workunit_per_cell_per_population( net.populations.contents.size() ); // will extend to include compartments, somehow LATER
	
	
	// Since the model is presented in its entirety, nodes need to perform domain decomposition themselves
	// Unfortunately, this means all nodes need to consider the existence of neuron Global ID's of all peers
	// LATER make a file format that enables fully distributed loading, through pre-processed domain decomposition (using e.g. METIS, or Scotch)
	
	// Mapping of neuron GIDs <-> ( nodes, work items, PointOnCellLocator's )
	
	std::map< Int, int > neuron_gid_to_node;	
	std::vector< std::map<Int, Int> > neuron_gid_per_cell_per_population( net.populations.contents.size() );
	
	// for local neurons only
	std::map< Int, work_t > neuron_gid_to_workitem;
	
	// similar to PointOnCellLocator
	// TODO refactor better
	struct CellLocator_PopInst {
		Int pop_seq;
		Int inst_seq;
		CellLocator_PopInst( Int _p, Int _i ){
			pop_seq = _p;
			inst_seq = _i;
		}
	};
	std::map< Int, CellLocator_PopInst > neuron_gid_to_popinst;
	
	// and helper functions, to not mess with the data structures directly (TODO refactor into object)
	auto GetLocalWorkItem_FromPopInst = [ &net, &local_workunit_per_cell_per_population  ]( Int pop_seq, Int cell_seq ){
		
		// sanity check
		if(!( 0 <= pop_seq && pop_seq < (Int)net.populations.contents.size() )) return (ptrdiff_t) -1;
		
		auto &hm = local_workunit_per_cell_per_population[pop_seq];
		
		if( !hm.count(cell_seq) ) return (ptrdiff_t) -1;
		
		return hm.at(cell_seq);
	};
	
	auto GetGlobalGid_FromPopInst = [ &net, &neuron_gid_per_cell_per_population  ]( Int pop_seq, Int cell_seq ){
		
		// sanity check
		if(!( 0 <= pop_seq && pop_seq < (Int)net.populations.contents.size() )) return (ptrdiff_t) -1;
		
		auto &hm = neuron_gid_per_cell_per_population[pop_seq];
		
		if( !hm.count(cell_seq) ) return (ptrdiff_t) -1;
		
		return hm.at(cell_seq);
	};
	
	auto GetRemoteNode_FromPopInst = [ &GetGlobalGid_FromPopInst, &neuron_gid_to_node ]( Int pop_seq, Int cell_seq ){
		// sanity check
		
		Int gid = GetGlobalGid_FromPopInst( pop_seq, cell_seq );
		if( gid < 0 ) return (int) ~0xABadD00d;
		
		if( !neuron_gid_to_node.count(gid) ){
			printf("Internal error: missing node for neuron gid %ld\n", gid);
			return (int) ~0xABadD00d;
		}
		
		return neuron_gid_to_node.at(gid);
	};
	
	// Maps neuron instance to either non-negative local work item, or negative ~(remote_node_id)
	auto WorkUnitOrNode = [ &GetLocalWorkItem_FromPopInst, &GetRemoteNode_FromPopInst ]( int pop, int cell_inst ){
		work_t ret = GetLocalWorkItem_FromPopInst( pop, cell_inst );
		// Say("pop %d %d = %llx", pop, cell_inst, (long long)ret);
		
		if( ret < 0 ){
			ret = ~GetRemoteNode_FromPopInst( pop, cell_inst );
		}
		// Say("popp %d %d = %llx", pop, cell_inst, (long long)ret);
		
		return ret;
	};
	
	#else
	
	// Auxiliary map to retrieve each work unit instance's position in the global tables
	// It's useful to remember that 'index' tabs should map contents for each work unit,
	// 	(and thus thould have the same length as associated data tables!)
	// so only a realized population -> work items mapping needs to be maintained
	
	std::vector< std::vector<size_t> >  workunit_per_cell_per_population( net.populations.contents.size() ); // will extend to include compartments, somehow LATER
	// GetLocalWorkItem_FromPopInst
	
	#endif
	
	printf("Creating populations...\n");
	
	auto InstantiateCellAsWorkitem = [ &config, &input_sources, &tabs ](
		const CellType &cell_type, const CellInternalSignature &sig,
		Int cell_gid, // for intra-cell randomization
		Int simulation_rng_seed,
		size_t &work_unit
	){
		(void) config; // just in case
		
		const auto &wig = sig.cell_wig;		
		
		//keep where the workitem starts
		work_unit = tabs.callbacks.size();
		
		// instantiate internal working sets
		size_t local_state_f32_index = tabs.global_initial_state.size();
		tabs.global_state_f32_index.push_back(local_state_f32_index);
		AppendToVector(tabs.global_initial_state, wig.state);
		
		size_t local_const_f32_index = tabs.global_constants.size();
		tabs.global_const_f32_index.push_back(local_const_f32_index);
		AppendToVector(tabs.global_constants, wig.constants);
		
		// and populate any scalar constants
		
		// the RNG seed for the cell
		ptrdiff_t Index_RngSeed = sig.common_in_cell.cell_rng_seed.Index_RngSeed;
		if( Index_RngSeed >= 0 ){
			
			// TODO use all bits
			// TODO use a more convenient/effective/powerful/etc algorithm to pass the simulation seed 
			// This one prevents a low seed from interfering with a low neuron gid, with no collision risk in XOR mixing
			
			auto ReverseBits = []( uint32_t x ) {
				x = (x & 0xFFFF0000) >> 16 | (x & 0x0000FFFF) << 16;
				x = (x & 0xFF00FF00) >>  8 | (x & 0x00FF00FF) <<  8;
				x = (x & 0xF0F0F0F0) >>  4 | (x & 0x0F0F0F0F) <<  4;
				x = (x & 0xCCCCCCCC) >>  2 | (x & 0x33333333) <<  2;
				x = (x & 0xAAAAAAAA) >>  1 | (x & 0x55555555) <<  1;
				return x;
			};
			
			// TODO use all bits in RNG seed
			uint32_t combined_seed = ReverseBits( (uint32_t) simulation_rng_seed ) ^ (uint32_t) cell_gid;
			
			tabs.global_constants[ local_const_f32_index + Index_RngSeed ] = EncodeI32ToF32( (int32_t) combined_seed );
		}
		
		
		// instantiate extension tables
		// also perhaps invert order of instantiantion/population and/or use heuristics(inputs, synapses, density...) to make the decision LATER
		auto AppendNewTables = [](auto &global_table_index, auto &global_table_arrays,  const size_t this_many){
			global_table_index.push_back( global_table_arrays.size() );
			for( size_t i = 0; i < this_many; i++){
				global_table_arrays.emplace_back();
			}
		};
		AppendNewTables(tabs.global_table_const_f32_index, tabs.global_tables_const_f32_arrays, wig.tables_const_f32.size());
		AppendNewTables(tabs.global_table_const_i64_index, tabs.global_tables_const_i64_arrays, wig.tables_const_i64.size());
		AppendNewTables(tabs.global_table_state_f32_index, tabs.global_tables_state_f32_arrays, wig.tables_state_f32.size());
		AppendNewTables(tabs.global_table_state_i64_index, tabs.global_tables_state_i64_arrays, wig.tables_state_i64.size());
		
		// and also populate the tables
		
		const auto off_cf32 = tabs.global_table_const_f32_index[work_unit];
		const auto off_sf32 = tabs.global_table_state_f32_index[work_unit];
		const auto off_ci64 = tabs.global_table_const_i64_index[work_unit];
		auto &tab_cf32 = tabs.global_tables_const_f32_arrays;
		auto &tab_sf32 = tabs.global_tables_state_f32_arrays;
		auto &tab_ci64 = tabs.global_tables_const_i64_arrays;
		
		if( cell_type.type == CellType::PHYSICAL ){
			// TODO profile and investigate instantiation performance, AppendToVector
			const auto &pig = sig.physical_cell;
			// tables for cell equation solver TODO revise when per-compartment analysis is done
			
			// adjacent compartment vectors for non-flattened, localized cable equation (like fwd or diagonal)
			for( size_t comp_seq = 0; comp_seq < pig.comp_implementations.size(); comp_seq++ ){
				auto &comp_impl = pig.comp_implementations[comp_seq];
				if( comp_impl.Index_AdjComp >= 0 ){
					RawTables::Table_I64 &AdjCompList = tab_ci64[ off_ci64 + comp_impl.Index_AdjComp ];
					AppendToVector( AdjCompList, pig.comp_definitions[comp_seq].adjacent_compartments );
				}
			}
			
			// remap vectors for compartment code deduplication
			if( pig.compartment_grouping == CellInternalSignature::CompartmentGrouping::GROUPED ){
				
				const auto &gp = pig.comp_group_impl;
				
				for( std::size_t comptype_seq = 0; comptype_seq < gp.distinct_compartment_types.size(); comptype_seq++ ){
					RawTables::Table_I64 &CompList = tab_ci64[ off_ci64 + gp.Index_CompList[comptype_seq] ];
					AppendToVector( CompList, gp.distinct_compartment_types[comptype_seq].toArray() );
				}
				
				RawTables::Table_I64 &Roff    = tab_ci64[off_ci64 + gp.Index_Roff   ]; AppendToVector( Roff    , gp.r_off    );
				RawTables::Table_I64 &Coff    = tab_ci64[off_ci64 + gp.Index_Coff   ]; AppendToVector( Coff    , gp.c_off    );
				RawTables::Table_I64 &Soff    = tab_ci64[off_ci64 + gp.Index_Soff   ]; AppendToVector( Soff    , gp.s_off    );
				RawTables::Table_I64 &CF32off = tab_ci64[off_ci64 + gp.Index_CF32off]; AppendToVector( CF32off , gp.cf32_off    );
				RawTables::Table_I64 &SF32off = tab_ci64[off_ci64 + gp.Index_SF32off]; AppendToVector( SF32off , gp.sf32_off    );
				RawTables::Table_I64 &CI64off = tab_ci64[off_ci64 + gp.Index_CI64off]; AppendToVector( CI64off , gp.ci64_off    );
				RawTables::Table_I64 &SI64off = tab_ci64[off_ci64 + gp.Index_SI64off]; AppendToVector( SI64off , gp.si64_off    );
				
			}
			
			// also populate tables of backward solvers (perhaps in a different loop LATER?)
			const auto &cabl_impl = pig.cable_solver_implementation;
			const auto &cabl_def = pig.cable_solver;
			if( cabl_def.type == SimulatorConfig::CABLE_FWD_EULER ){
				// no helper arrays
			}
			else if( cabl_def.type == SimulatorConfig::CABLE_BWD_EULER ){
				RawTables::Table_I64 &Order  = tab_ci64[off_ci64 + cabl_impl.Index_BwdEuler_OrderList];
				RawTables::Table_I64 &Parent = tab_ci64[off_ci64 + cabl_impl.Index_BwdEuler_ParentList];
				RawTables::Table_F32 &InvRCD = tab_cf32[off_cf32 + cabl_impl.Index_BwdEuler_InvRCDiagonal];
				RawTables::Table_F32 &WorkD  = tab_sf32[off_sf32 + cabl_impl.Index_BwdEuler_WorkDiagonal];
				
				AppendToVector(Order, cabl_def.BwdEuler_OrderList);
				AppendToVector(Parent, cabl_def.BwdEuler_ParentList);
				AppendToVector(InvRCD, cabl_def.BwdEuler_InvRCDiagonal);
				WorkD.resize(Order.size(), NAN); // don't care about the content, but must initialize it
			}
			else{
				printf("Unknown cable solver %d for %s\n", cabl_def.type, sig.name.c_str());
				return false;
			}
			
		}
		
		if( cell_type.type == CellType::ARTIFICIAL ){
			const auto &cell = cell_type.artificial;
			const auto &aig = sig.artificial_cell;
			
			if( cell.type == ArtificialCell::SPIKE_SOURCE ){
				const InputSource &input = input_sources.get(cell.spike_source_seq);
				const auto &inpimp = aig.inpimpl;
				
				if( input.type == InputSource::SPIKE_LIST ){
					// the spike list needs its own array
					// append spike list to the common vector, also a sentinel to avoid checking the indices
					RawTables::Table_F32 &times = tab_cf32[off_cf32 + inpimp.Table_SpikeListTimes];
					for( auto spike : input.spikes ) times.push_back( spike.time_of_occurrence );
					times.push_back( FLT_MAX ); // sentinel value
				}
				
			}
		}
		
		// instantiate iteration callback
		tabs.callbacks.push_back(sig.callback);
		
		return true;
	};
	
	#ifdef USE_MPI
	// For domain decompoosition, get the footprint of the populations
	// like total amount of cells, etc.
	
	int total_neurons = 0; // LATER replace with referrable parts of the model, or sth
	for( Int pop_seq = 0; pop_seq < (Int)net.populations.contents.size() ; pop_seq++ ){
		const Network::Population &pop = net.populations.contents[pop_seq];
		total_neurons += (int) pop.instances.size();
	}
	
	Say("Total neurons: %d", total_neurons);
	
	
	// TODO this is where splitting happens
	// TODO simplest domain decomposition: by neuron GID
	// TODO domain decomposition shoould be done as a preliminary step,
	// 	so each node bothers only with its own part (and adjacent)
	// 	instead of trying to find out what is own, which work items and nodes are adjacent
	struct NodeMapper{
		// evened 
		// LATER more strategies
		int total_nodes;
		int total_items;
		
		NodeMapper(int _no, int _ne){
			total_nodes = _no;
			total_items = _ne;
		}
		
		int GetNodeFor(int item_gid ){
			
			// sanity check
			if(!( 0 <= item_gid && item_gid < total_items )) return -1;
			
			// distribute the items, as even as possible
			int work_items_evenly  = total_items / total_nodes;
			int work_items_residue = total_items % total_nodes;
			// each node gets a contiguous slice of GID's, for now
			// with lengths [ evenly + 1, evenly + 1 ... evenly, evenly ]
			
			int work_items_from_resi_nodes = (work_items_evenly + 1) * work_items_residue;
			if( item_gid < work_items_from_resi_nodes ){
				// in evenly + 1 part of nodes
				return ( item_gid ) / ( work_items_evenly + 1 );
			}
			else{
				// in part of nodes without the residue
				int off_nonres = item_gid - work_items_from_resi_nodes;
				return work_items_residue + ( off_nonres ) / ( work_items_evenly );
			}
		}
	};
	
	NodeMapper to_node( my_mpi.world_size, total_neurons );
	#endif
	
	timeval time_pops_start, time_pops_end;
	gettimeofday(&time_pops_start, NULL);
	
	int current_neuron_gid = 0;
	for( Int pop_seq = 0; pop_seq < (Int)net.populations.contents.size() ; pop_seq++ ){
		const Network::Population &pop = net.populations.contents[pop_seq];
		
		const CellType &cell_type = model.cell_types.get(pop.component_cell);
		const auto &sig = cell_sigs[pop.component_cell];
		
		// TODO pre-allocate since the values will be cloned in a predictable pattern
		for( Int inst_seq = 0; inst_seq < (Int)pop.instances.size(); inst_seq++ ){
			
			bool instantiate_this = true;
			
			#ifdef USE_MPI
			
			int on_node = to_node.GetNodeFor( current_neuron_gid );
			
			neuron_gid_per_cell_per_population[pop_seq][inst_seq] = current_neuron_gid;
			neuron_gid_to_node[current_neuron_gid] = on_node;
			
			if( on_node == my_mpi.rank ) instantiate_this = true;
			else instantiate_this = false;
			
			#endif
			
			if( instantiate_this ){
				size_t work_unit = SIZE_MAX;
				// if(config.debug){
				// printf("Instantiating cell %d...\n", current_neuron_gid);
				// }
				if( !InstantiateCellAsWorkitem( cell_type, sig, current_neuron_gid, simulation_random_seed, work_unit ) ) return false;
				
				#ifdef USE_MPI
				
				if(config.debug_netcode) Say("Instantiate %d %d -> %d", pop_seq, inst_seq, current_neuron_gid );
				
				local_workunit_per_cell_per_population[pop_seq][inst_seq] = work_unit;
				neuron_gid_to_workitem[current_neuron_gid] = work_unit;
				neuron_gid_to_popinst.insert( std::make_pair( current_neuron_gid, CellLocator_PopInst(pop_seq, inst_seq) ) );
				
				#else
				
				workunit_per_cell_per_population[pop_seq].push_back(work_unit);
				
				#endif
			}
			
			current_neuron_gid++;
		}
		
	}
	gettimeofday(&time_pops_end, NULL);
	printf("Created populations in %.4lf sec.\n",TimevalDeltaSec(time_pops_start, time_pops_end));
	
	// Add some extra misc-purpose tables
	tabs.global_const_tabref = tabs.global_tables_const_f32_arrays.size();
	tabs.                           global_tables_const_f32_arrays.emplace_back();
	tabs.global_state_tabref = tabs.global_tables_state_f32_arrays.size();
	tabs.                           global_tables_state_f32_arrays.emplace_back();
	
	// tabs.logger_table_const_i64_to_state_f32_index = tabs.global_tables_const_i64_arrays.size();
	// tabs.global_tables_const_f32_arrays.emplace_back();
	
	// Now add the attachments, as table entries
	
	// instantiation of synapse internals is also used in firing-synapse inputs
	auto AppendSyncompInternals = [ &DescribeLems_AppendTableEntry ](
		const SynapticComponent &syn, Int id_id, size_t work_unit, const CellInternalSignature::SynapticComponentImplementation &synimpl,
		RawTables &tabs
	){
		
		// no integers used for now
		const auto off_cf32 = tabs.global_table_const_f32_index[work_unit];
		const auto off_sf32 = tabs.global_table_state_f32_index[work_unit];
		auto &tab_cf32 = tabs.global_tables_const_f32_arrays;
		auto &tab_sf32 = tabs.global_tables_state_f32_arrays;
		
		// TODO eliminate id_id from here
		if(id_id < 0){
			SynapticComponent::Type core_id = SynapticComponent::Type(id_id + SynapticComponent::Type::MAX);
			switch(core_id){
			
			case SynapticComponent::Type::EXP :{
				
				RawTables::Table_F32 &Gbase = tab_cf32.at(off_cf32 + synimpl.Table_Gbase);
				RawTables::Table_F32 &Erev  = tab_cf32.at(off_cf32 + synimpl.Table_Erev );
				RawTables::Table_F32 &Tau   = tab_cf32.at(off_cf32 + synimpl.Table_Tau );
				
				RawTables::Table_F32 &Grel  = tab_sf32.at(off_sf32 + synimpl.Table_Grel );							
			
				Gbase.push_back(syn.exp.gbase);
				Erev .push_back(syn.exp.erev);
				Tau  .push_back(syn.exp.tauDecay);
				
				Grel.push_back(0); // LATER initialize the synapses somehow
				
				//could re-use globals LATER
				break;
			}
			case SynapticComponent::Type::GAP :{
				
				RawTables::Table_F32 &Gsyn = tab_cf32.at(off_cf32 + synimpl.Table_Gbase);
		
				Gsyn.push_back(syn.gap.conductance);
				
				break;
			}
			default:
				printf("internal error: populate unknown syncomp core_id %d\n", core_id);
				return false;
			}
		}
		else{
			if(syn.type == SynapticComponent::Type::BLOCKING_PLASTIC ){
				if( syn.blopla.block_mechanism.type != SynapticComponent::BlockingPlasticSynapse::BlockMechanism::NONE ){
					DescribeLems_AppendTableEntry( work_unit, syn.blopla.block_mechanism.component, synimpl.block_component );
				}
					
				if( syn.blopla.plasticity_mechanism.type != SynapticComponent::BlockingPlasticSynapse::PlasticityMechanism::NONE ){	
					DescribeLems_AppendTableEntry( work_unit, syn.blopla.plasticity_mechanism.component, synimpl.plasticity_component );
				}
				
				DescribeLems_AppendTableEntry( work_unit, syn.component, synimpl.synapse_component );
			}
			else if( syn.component.ok() ){
				DescribeLems_AppendTableEntry( work_unit, syn.component, synimpl.synapse_component );
			}
			else{
				printf("internal error: populate unknown syncomp id %ld\n", id_id);
				return false;
			}
		}
		return true;
	};
	
	// TODO get direct access to cell type or comp.implementation, could be wasteful to access cell_types repeatedly when cell_type is already known; or profile at least
	// TODO also, this smells like the properties could be refactored into compact objects
	auto GetCompartmentInputImplementations = [ &cell_types ]( const CellInternalSignature &sig, Int celltype_seq, Int seg_seq, Real fractionAlong ){
		const CellType &cell_type = cell_types.get(celltype_seq);
		if( cell_type.type == CellType::PHYSICAL ){
			auto comp_seq = sig.physical_cell.compartment_discretization.GetCompartmentForSegmentLocation(seg_seq, fractionAlong);
			return sig.physical_cell.comp_implementations.at(comp_seq).input;
		}
		else{
			return sig.artificial_cell.input;
		}
		
	};
	
	auto GetCompartmentSynapseImplementations = [ &cell_types, &cell_sigs, &net ]( const PointOnCellLocator &loc ){
		
		const auto &pop = net.populations.get(loc.population);
		// TODO for split cell
		auto celltype_seq = pop.component_cell;
		const CellInternalSignature &sig = cell_sigs.at(celltype_seq);
		const CellType &cell_type = cell_types.get(celltype_seq);
		
		if( cell_type.type == CellType::PHYSICAL ){
			auto comp_seq = sig.physical_cell.compartment_discretization.GetCompartmentForSegmentLocation(loc.segment, loc.fractionAlong);
			return sig.physical_cell.comp_implementations.at(comp_seq).synapse;
		}
		else{
			return sig.artificial_cell.synapse;
		}
		
	};
	auto GetCompartmentSpikerImplementation = [ &cell_types ]( const CellInternalSignature &sig, Int celltype_seq, Int seg_seq, Real fractionAlong ){
		const CellType &cell_type = cell_types.get(celltype_seq);
		if( cell_type.type == CellType::PHYSICAL ){
			auto comp_seq = sig.physical_cell.compartment_discretization.GetCompartmentForSegmentLocation(seg_seq, fractionAlong);
			return sig.physical_cell.comp_implementations.at(comp_seq).spiker;
		}
		else{
			return sig.artificial_cell.spiker;
		}
		
	};
	
	auto GetCompartmentVoltageStatevarIndex = [ &cell_types ]( const CellInternalSignature &sig, Int celltype_seq, Int seg_seq, Real fractionAlong ){
		const CellType &cell_type = cell_types.get(celltype_seq);
		if( cell_type.type == CellType::PHYSICAL ){
			auto comp_seq = sig.physical_cell.compartment_discretization.GetCompartmentForSegmentLocation(seg_seq, fractionAlong);
			return (ptrdiff_t) sig.physical_cell.GetVoltageStatevarIndex( comp_seq );
		}
		else{
			return sig.artificial_cell.Index_Statevar_Voltage;
		}
	};
	
	#ifdef USE_MPI
	// Variants of GetXxxImplementation, including ???
	auto GetCompartmentSpikerImplementation_Global = [ &net, &cell_sigs, &tabs, &GetLocalWorkItem_FromPopInst, &GetCompartmentSpikerImplementation ]( const PointOnCellLocator &loc, size_t &spiker_table_idx ){
		
		const auto &pop = net.populations.get(loc.population);
		// TODO for split cell
		auto celltype_seq = pop.component_cell;
		const CellInternalSignature &sig = cell_sigs.at(celltype_seq);
		
		ptrdiff_t work_unit = GetLocalWorkItem_FromPopInst( loc.population, loc.cell_instance );
		assert( work_unit >= 0 ); // LATER make fallible ?
		
		ptrdiff_t local_offset = GetCompartmentSpikerImplementation( sig, celltype_seq, loc.segment, loc.fractionAlong ).Table_SpikeRecipients;
		assert( local_offset >= 0 );
		
		spiker_table_idx = tabs.global_table_const_i64_index[work_unit] + local_offset;
		
		return true;
	};
	auto GetCompartmentVoltageStatevarIndex_Global = [ &net, &cell_sigs, &tabs, &GetLocalWorkItem_FromPopInst, &GetCompartmentVoltageStatevarIndex ]( const PointOnCellLocator &loc ){
		const auto &pop = net.populations.get(loc.population);
		// TODO for split cell
		auto celltype_seq = pop.component_cell;
		const CellInternalSignature &sig = cell_sigs.at(celltype_seq);
		
		ptrdiff_t work_unit = GetLocalWorkItem_FromPopInst( loc.population, loc.cell_instance );
		assert( work_unit >= 0 ); // LATER make fallible ?
		
		ptrdiff_t local_offset = GetCompartmentVoltageStatevarIndex( sig, celltype_seq, loc.segment, loc.fractionAlong );
		assert( local_offset >= 0 );
		
		size_t global_idx_V_peer = tabs.global_state_f32_index[work_unit] + local_offset;
		return global_idx_V_peer;
	};
	#endif
	// create connectivity after cells are instantiated, for simplicity:
	
	// also populate the inputs
	printf("Creating inputs...\n");
	
	timeval time_inps_start, time_inps_end;
	gettimeofday(&time_inps_start, NULL);
	
	for(size_t inp_seq = 0; inp_seq < net.inputs.size(); inp_seq++){
		const auto &inp = net.inputs[inp_seq];
		// printf("input cell %ld seq %ld\n", (Int)inp.cell_instance, (Int)inp.segment );
		const auto &source = input_sources.get(inp.component_type);
		const auto &pop = net.populations.get(inp.population);
		
		// get Cell type/Compartment instance
		const auto &sig = cell_sigs[pop.component_cell];
		
		//get work unit from type/population/instance
		#ifdef USE_MPI
		
		ptrdiff_t work_unit = GetLocalWorkItem_FromPopInst( inp.population, inp.cell_instance );
		if( work_unit < 0 ) continue;
		
		#else
		
		size_t work_unit = workunit_per_cell_per_population[inp.population][inp.cell_instance];
		// branch for whole cell or compartment LATER
		#endif
		
		const auto &inpimps = GetCompartmentInputImplementations(sig, pop.component_cell, inp.segment, inp.fractionAlong);
		
		// get Input type
		Int id_id = GetInputIdId( inp.component_type );
		
		if( !inpimps.count(id_id) ){
			printf("Internal error: No input implementation for input type %ld\n", id_id);
			return false;
		}
		const CellInternalSignature::InputImplementation &inpimp = inpimps.at(id_id);
		
		// consider the common tables such as weight
		const auto off_cf32 = tabs.global_table_const_f32_index[work_unit];
		// const auto off_sf32 = tabs.global_table_state_f32_index[work_unit];
		const auto off_si64 = tabs.global_table_state_i64_index[work_unit];
		auto &tab_cf32 = tabs.global_tables_const_f32_arrays;
		// auto &tab_sf32 = tabs.global_tables_state_f32_arrays;
		auto &tab_si64 = tabs.global_tables_state_i64_arrays;
		
		// TODO something about the LEMS weight property, or ignore it if LEMS does so
		float weight = inp.weight;
		if( !std::isfinite(weight) ) weight = 1;
		tab_cf32[ off_cf32 + inpimp.Table_Weight ].push_back(weight);
		
		// helpers
		auto PopulateSpikeList = [ &tab_cf32, &off_cf32, &tab_si64, &off_si64 ]( const auto &spike_list, auto &inpimp ){
			// append spike list to the common vector, also a sentinel to avoid checking the indices
			// also append a start and an initial position
			RawTables::Table_F32 &times = tab_cf32[off_cf32 + inpimp.Table_SpikeListTimes];
			//RawTables::Table_F32 &starts = tab_ci64[off_ci64 + inpimp.Table_SpikeListStarts];
			RawTables::Table_I64 &positions = tab_si64[off_si64 + inpimp.Table_SpikeListPos];
			
			positions.push_back( times.size() );
			
			for( auto spike : spike_list ) times.push_back( spike.time_of_occurrence );
			times.push_back( FLT_MAX ); // sentinel value
			
		};
		
		// could re-use globals LATER
		if(id_id < 0){
			InputSource::Type core_id = InputSource::Type(id_id + InputSource::Type::MAX);
			
			switch(core_id){
				
			case InputSource::Type::PULSE :{
				
				RawTables::Table_F32 &Imax = tab_cf32[off_cf32 + inpimp.Table_Imax];
				RawTables::Table_F32 &start = tab_cf32[off_cf32 + inpimp.Table_Delay];
				RawTables::Table_F32 &duration = tab_cf32[off_cf32 + inpimp.Table_Duration];
				
				Imax.push_back(source.amplitude);
				start.push_back(source.delay);
				duration.push_back(source.duration);
				
				break;
			}
			case InputSource::Type::SPIKE_LIST :{
				
				PopulateSpikeList( source.spikes, inpimp );
				
				break;
			}
			default:
				// internal error
				printf("populate: Unknown input core_id %d\n", core_id);
				return false;
			}
		}
		else{
			if(
				source.type == InputSource::Type::TIMED_SYNAPTIC
				|| source.type == InputSource::Type::POISSON_SYNAPSE
				|| source.type == InputSource::Type::POISSON_SYNAPSE_TRANSIENT
			){
				
				if( source.type == InputSource::Type::TIMED_SYNAPTIC ){
					PopulateSpikeList( source.spikes, inpimp );
				}
				else if(
					source.type == InputSource::Type::POISSON_SYNAPSE
					|| source.type == InputSource::Type::POISSON_SYNAPSE_TRANSIENT
				){
					DescribeLems_AppendTableEntry( work_unit, source.component, inpimp.component );
				}
				else{
					printf("internal error: input component %ld append for what sort of firing synapse input? \n", id_id);
					return false;
				}
				
				// and append comp. signature to tables of same type
				const auto &syn = synaptic_components.get(source.synapse);
				if( !AppendSyncompInternals( syn, GetSynapseIdId(source.synapse), work_unit, inpimp.synimpl, tabs ) ) return false;
				
			}
			else if( source.component.ok() ){
				DescribeLems_AppendTableEntry( work_unit, source.component, inpimp.component );
			}
			else{
				printf("internal error: populate unknown input id %ld\n", id_id);
				return false;
			}
		}
		
	}
	
	gettimeofday(&time_inps_end, NULL);
	printf("Created inputs in %.4lf sec.\n",TimevalDeltaSec(time_inps_start, time_inps_end));
	
	// also populate the synapses
	// place the append syncomp lambda somewhere here LATER
	printf("Creating synapses...\n");
	
	timeval time_syns_start, time_syns_end;
	gettimeofday(&time_syns_start, NULL);
	
	for(size_t proj_seq = 0; proj_seq < net.projections.contents.size(); proj_seq++){
		
		// printf("Projection %zd of %zd \n", proj_seq, net.projections.contents.size() );
		
		const auto &proj = net.projections.contents.at(proj_seq);
		const auto &prepop = net.populations.get(proj.presynapticPopulation);
		const auto &postpop = net.populations.get(proj.postsynapticPopulation);
		
		// get Cell type/Compartment instance. TODO remove.
		const auto &presig = cell_sigs[prepop.component_cell];
		const auto &postsig = cell_sigs[postpop.component_cell];
		
		auto AppendSynapticComponentEntries = [
			&model, &prepop, &AppendSyncompInternals, &GetSynapseIdId,
			&GetCompartmentSynapseImplementations, &GetCompartmentSpikerImplementation, &GetCompartmentVoltageStatevarIndex
			#ifdef USE_MPI
			, &AppendRemoteDependency_Vpeer, &AppendRemoteDependency_Spike
			#endif
		](
			const SynapticComponent &syn, Int syncomp_seq, const Network::Projection::Connection &conn,
			const PointOnCellLocator &mine_loc,
			const PointOnCellLocator &peer_loc,
			work_t work_unit, const CellInternalSignature &sig, Int mine_cell_type_seq,
			work_t peer_work_unit,const CellInternalSignature &peer_sig, Int peer_cell_type_seq,
			RawTables &tabs
		){
			// TODO perhaps early exit if local work item is remote
			// branch for whole cell or compartment LATER
			
			Int id_id = GetSynapseIdId( syncomp_seq );
			
			const bool needs_spike = syn.HasSpikeIn(model.component_types);
			const bool needs_Vpeer = syn.HasVpeer(model.component_types);
			
			// get weight value early on, just in case remote may need it ...?
			Real weight = conn.weight;
			if( !std::isfinite(weight) ) weight = 1;
			
			// get the underlying mechanism
			// TODO might not exist, if this node doesn't work with this cell type
			const auto &synimps = GetCompartmentSynapseImplementations( mine_loc );
			if(!synimps.count(id_id)){
				printf("Internal error: No impl signature for type %ld\n", id_id);
				printf("Synimps: " );
				for( auto keyval : synimps ) printf("%ld ", keyval.first );
				printf("\n" );
				return false;
			}
			const CellInternalSignature::SynapticComponentImplementation &synimpl = synimps.at(id_id);
			
			// also add weight
			auto AddWeight = [ &tabs ]( work_t work_unit, const auto &synimpl, Real weight ){
				if( work_unit < 0 ){
					return true; // not on this node
				}
				const auto off_cf32 = tabs.global_table_const_f32_index[work_unit];
				auto      &tab_cf32 = tabs.global_tables_const_f32_arrays;
				
				tab_cf32[ off_cf32 + synimpl.Table_Weight ].push_back(weight);
				return true;
			};
			bool uses_weight = true; // might be elided LATER
			if( uses_weight ) AddWeight( work_unit, synimpl, weight );
			
			if( needs_Vpeer ){
				
				auto AddGap = [ 
					&GetCompartmentVoltageStatevarIndex
					#ifdef USE_MPI
					, &AppendRemoteDependency_Vpeer
					#endif
				](
					const SynapticComponent &syn,
					work_t work_unit, const CellInternalSignature::SynapticComponentImplementation &synimpl, //const CellInternalSignature &sig, Int comp_seq,
					work_t peer_work_unit,const CellInternalSignature &peer_sig, Int peer_cell_type_seq,
					const PointOnCellLocator &peer_loc,					
					auto &tabs
				){
					
					if( work_unit < 0 ){
						return true; // not mine, let the need to send emerge
					}
					
					
					const auto off_ci64 = tabs.global_table_const_i64_index[work_unit];
					auto &tab_ci64 = tabs.global_tables_const_i64_arrays;
					
					auto glob_tab_Vpeer = off_ci64 + synimpl.Table_Vpeer;
					
					RawTables::Table_I64 &Vpeer = tab_ci64.at(glob_tab_Vpeer);
					
					// if it exists on this node
					#if USE_MPI
					if( peer_work_unit < 0 ){
						
						int node_peer = ~(peer_work_unit);
						// get the value from remote peer
						if( !AppendRemoteDependency_Vpeer( peer_loc, node_peer, glob_tab_Vpeer ) ) return false;
						return true;
					}
					#endif
					// otherwise it's local
					ptrdiff_t local_idx_V_peer = GetCompartmentVoltageStatevarIndex( peer_sig, peer_cell_type_seq, peer_loc.segment, peer_loc.fractionAlong );
						
					if( local_idx_V_peer < 0 ){
						printf("internal error: gap junction realization: Cell type %ld has no Vpeer\n",peer_cell_type_seq );
						return false;
					}
					
					// get reference to where peer's voltage is located
					ptrdiff_t global_idx_V_peer = tabs.global_state_f32_index[peer_work_unit] + local_idx_V_peer;
					auto global_tabentry = GetEncodedTableEntryId( tabs.global_state_tabref, global_idx_V_peer);
					
					Vpeer.push_back(global_tabentry);
					
					return true;
				};
				
				if( !AddGap(syn, work_unit, synimpl, peer_work_unit, peer_sig, peer_cell_type_seq, peer_loc, tabs) ) return false;
			}
			
			if( needs_spike ){
				// LATER validate conditions one more time:
				// pre-synaptic must have a spike output port
				// post-synaptic must have a spike input
				
				auto AddDelay = [ &tabs ]( work_t work_unit, const auto &synimpl, Real delay ){
					if( work_unit < 0 ){
						return true; // not on this node
					}
					const auto off_cf32 = tabs.global_table_const_f32_index[work_unit];
					auto      &tab_cf32 = tabs.global_tables_const_f32_arrays;
					const auto off_sf32 = tabs.global_table_state_f32_index[work_unit];
					auto      &tab_sf32 = tabs.global_tables_state_f32_arrays;
					
					tab_cf32[ off_cf32 + synimpl.Table_Delay ].push_back(delay);
					tab_sf32[ off_sf32 + synimpl.Table_NextSpike ].push_back(-INFINITY);
					
					return true;
				};
				
				// delay is common for all 
				bool uses_delay = true;
				if( uses_delay ){
					
					Real delay = conn.delay;
					if( !std::isfinite(conn.delay) ) delay = 0;
					
					AddDelay( work_unit, synimpl, delay );
				}
				
				
				auto AddChemPrePost = [
					&prepop, &GetCompartmentSpikerImplementation
					#ifdef USE_MPI
					, &AppendRemoteDependency_Spike
					#endif
				](
					const SynapticComponent &syn,
					work_t post_work_unit, const CellInternalSignature::SynapticComponentImplementation &post_synimpl,
					work_t pre_work_unit,const CellInternalSignature &pre_sig, Int pre_cell_type_seq,
					const PointOnCellLocator &pre_loc,
					auto &tabs
				){
					
					if( post_work_unit < 0 ){
						return true; // not mine, let it be resolved on spike-receiver demands later on
					}
					
					// and now add the entries to the post syn table
					auto AddPost = [ ]( auto &tabs, work_t work_unit, const auto &synimpl){
						
						const auto off_si64 = tabs.global_table_state_i64_index[work_unit];
						auto &tab_si64 = tabs.global_tables_state_i64_arrays;
						
						RawTables::Table_I64 &Trig  = tab_si64.at(off_si64 + synimpl.Table_Trig );
						Trig.push_back(0); // perhaps compress trigger table LATER
						
						return true;
					};
					
					//post has trig buf, gives idx to sender
					
					// and where the spike target is located in each post-synaptic work item
						
					//TODO change for split cell? just find the compartment responsible
					
					long long global_idx_T_dest_table = tabs.global_table_state_i64_index[post_work_unit] + post_synimpl.Table_Trig; 
					// printf("yyyyyy %lld %lld %lld\n\n\n", post_work_unit, tabs.global_table_state_i64_index[post_work_unit], global_idx_T_dest_table );
					long long entry_idx_T_dest = tabs.global_tables_state_i64_arrays[global_idx_T_dest_table].size(); // TODO change to reflect when handling mask, perhaps encapsulate
					auto packed_id = GetEncodedTableEntryId( global_idx_T_dest_table, entry_idx_T_dest );
					
					if( !AddPost( tabs, post_work_unit, post_synimpl) ) return false;
					
					#ifdef USE_MPI
					if( pre_work_unit < 0 ){
						
						// needs to ask for spike input from pre, and keep the map from packed received input to buf
						int node_pre = ~(pre_work_unit); // TODO refactor
						if( !AppendRemoteDependency_Spike( pre_loc, node_pre, packed_id ) ) return false;
						return true;
					}
					#endif
					// local pre, needs idx from sender
						
					// XXX do this for emergent remote ones too
					
					// Spike output is required to exist on pre compartment
					const auto &preimp = GetCompartmentSpikerImplementation( pre_sig, pre_cell_type_seq, pre_loc.segment, pre_loc.fractionAlong );
					if( preimp.Table_SpikeRecipients < 0 ){
						printf("Internal error: No spike send for celltype %ld seg %ld fractionAlong %g %zd\n", prepop.component_cell, pre_loc.segment, pre_loc.fractionAlong, preimp.Table_SpikeRecipients);
						return false;
					}
					
					RawTables::Table_I64 &Spike_recipients = tabs.global_tables_const_i64_arrays.at(tabs.global_table_const_i64_index.at(pre_work_unit) + preimp.Table_SpikeRecipients);
					
					Spike_recipients.push_back(packed_id);
					
					return true;
				};
				
				if( !AddChemPrePost(syn, work_unit, synimpl, peer_work_unit, peer_sig, peer_cell_type_seq, peer_loc, tabs) ) return false;
			};
			
			if( work_unit >= 0 ){
				// and join them, implementing internal mech. just once
				if( !AppendSyncompInternals( syn, id_id, work_unit, synimpl, tabs ) ) return false;
			}
			
			
			return true; // yay!
		};
		
		for(size_t conn_seq = 0; conn_seq < proj.connections.contents.size(); conn_seq++){
			const auto &conn = proj.connections.contents[conn_seq];
			
			// printf("Connection %zd of %zd \n", conn_seq, proj.connections.contents.size() ); fflush(stdout);
			
			const PointOnCellLocator pre_loc  = { proj.presynapticPopulation , conn.preCell , conn.preSegment , conn.preFractionAlong  };
			const PointOnCellLocator post_loc = { proj.postsynapticPopulation, conn.postCell, conn.postSegment, conn.postFractionAlong };
			
			//get work unit from type/population/instance. TODO replace with something better, using work_unit only when applicable. (eg split cell)
			
			#ifdef USE_MPI
			
			work_t work_unit_pre  = WorkUnitOrNode( proj.presynapticPopulation , conn.preCell  );
			work_t work_unit_post = WorkUnitOrNode( proj.postsynapticPopulation, conn.postCell );
			
			#else 
			
			work_t work_unit_pre = workunit_per_cell_per_population[proj.presynapticPopulation][conn.preCell];
			work_t work_unit_post = workunit_per_cell_per_population[proj.postsynapticPopulation][conn.postCell];
			
			#endif
			
			// printf("connnn %lx %ld\n", work_unit_pre, work_unit_post); fflush(stdout);
			
			// now populate the synaptic components in the raw tables
			if(conn.type == Network::Projection::Connection::SPIKING){
				// TODO validate conditions:
				// pre-synaptic must have a spike output port
				// post-synaptic must have a spike input
				
				const SynapticComponent &syn = synaptic_components.get(conn.synapse);
				
				if( !AppendSynapticComponentEntries(
					syn, conn.synapse, conn,
					post_loc, pre_loc,
					work_unit_post, postsig, postpop.component_cell,
					work_unit_pre , presig , prepop.component_cell ,
					tabs
				) ) return false;
			}
			else if(conn.type == Network::Projection::Connection::ELECTRICAL){
				// same behaviour for pre-and post-synaptic
				// Vpeer is required, and it always exists for physical compartments
				
				const SynapticComponent &syn = synaptic_components.get(conn.synapse);
				
				if( !AppendSynapticComponentEntries(
					syn, conn.synapse, conn,
					post_loc, pre_loc,
					work_unit_post, postsig, postpop.component_cell,
					work_unit_pre , presig , prepop.component_cell ,
					tabs
				) ) return false;
				if( !AppendSynapticComponentEntries(
					syn, conn.synapse, conn,
					pre_loc, post_loc,
					work_unit_pre , presig , prepop.component_cell ,
					work_unit_post, postsig, postpop.component_cell,
					tabs
				) ) return false;	
			}
			else if(conn.type == Network::Projection::Connection::CONTINUOUS){
				// anything goes, really
				const SynapticComponent &syn_pre  = synaptic_components.get(conn.continuous.preComponent );
				const SynapticComponent &syn_post = synaptic_components.get(conn.continuous.postComponent);
				
				if( !AppendSynapticComponentEntries(
					syn_post, conn.continuous.postComponent, conn,
					post_loc, pre_loc,
					work_unit_post, postsig, postpop.component_cell,
					work_unit_pre , presig , prepop.component_cell ,
					tabs
				) ) return false;
				if( !AppendSynapticComponentEntries(
					syn_pre , conn.continuous.preComponent , conn,
					pre_loc, post_loc,
					work_unit_pre , presig , prepop.component_cell ,
					work_unit_post, postsig, postpop.component_cell,
					tabs
				) ) return false;
			}
			else{
				printf("internal error: populate unknown synapse type projection %zd instance %zd\n", proj_seq, conn_seq);
				return false;
			}
			// printf("conndone %zd %zd\n", proj_seq, conn_seq); fflush(stdout);
		}
		
	}
	
	gettimeofday(&time_syns_end, NULL);
	printf("Created synapses in %.4lf sec.\n",TimevalDeltaSec(time_syns_start, time_syns_end));
	
	printf("Creating data outputs...\n");
	
	// add the loggers
	
	// TODO explicit capture
	auto Implement_LoggerColumn = [ & ]( const Int daw_seq, const Int col_seq, const std::string &output_filepath,const auto &path, EngineConfig::TrajectoryLogger::LogColumn &column ){
		
		if( path.RefersToCell() ){
			
			const auto &pop = net.populations.get(path.population);
			const auto &cell_type = cell_types.get(pop.component_cell);
			
			// get work item for this neuron here, or perhaps further on when compartments are work items
			
			#ifdef USE_MPI
			
			work_t work_unit = WorkUnitOrNode( path.population, path.cell_instance );
			
			if( work_unit < 0 ){
				#ifdef MPI
				assert( my_mpi.rank == 0 );
				#endif
				
				int remote_node = ~(work_unit);
				column.on_node = remote_node;
				if( !AppendRemoteDependency_DataWriter( {daw_seq, col_seq}, remote_node ) ) return false;
				
				return true;
			}
			
			#else
			
			work_t work_unit = workunit_per_cell_per_population[path.population][path.cell_instance];
			
			#endif
			
			// get sig for said neuron, to properly manipulate the work item
			const auto &sig = cell_sigs[pop.component_cell];
			
			// commonly used paths (evidence mounts to use parent class form commonish things like inputs and synapses)
			
			// otherwise cell-type-specific paths
			
			auto MustBePhysicalCell = [ &col_seq, &output_filepath ]( const CellType &cell_type ){
				if( cell_type.type != CellType::PHYSICAL ){
					printf("internal error: column %ld for data writer %s has channel path on non-physical cell\n", col_seq, output_filepath.c_str());
					return false;
				}
				return true;
			};
			auto MustBeArtificialCell = [ &col_seq, &output_filepath ]( const CellType &cell_type ){
				if( cell_type.type != CellType::ARTIFICIAL ){
					printf("internal error: column %ld for data writer %s has cell path on non-artificial cell\n", col_seq, output_filepath.c_str());
					return false;
				}
				return true;
			};
			
			// NB currently, all cell-facing references have a valid reference to segment. (That is 0 for artificial cells.)
			// But this could not always hold LATER...
			// This is moved in this scope tentatively, to avoid copy pasting deeper in
			int32_t comp_seq = -1;
			if( cell_type.type == CellType::PHYSICAL ){
				if( path.segment_seq >= 0 ){
					comp_seq = GetCompSeqForCell(cell_type, sig, path.segment_seq, 0.5); // NB: LEMS refs lack fractionAlong for now
				}
				else{
					assert(false); // should not happen as of now, guard for when this case happens LATER
				}
			}
			
			switch(path.type){
			case Simulation::LemsQuantityPath::Type::SEGMENT: {
				
				if( !MustBePhysicalCell( cell_type ) ) return false;
				
				const PhysicalCell &cell = cell_type.physical;
				// const Morphology &morph = morphologies.get(cell.morphology);
				const BiophysicalProperties &bioph = biophysics.at(cell.biophysicalProperties);
				
				auto &pig = sig.physical_cell;
				
				// NB: If the path refers to a segment, 
				
				switch(path.segment.type){
					
					case Simulation::LemsQuantityPath::SegmentPath::Type::VOLTAGE: {
						// printf("daw volt cell %ld seg %ld comp %d\n", (Int)path.cell_instance, (Int)path.segment_seq, (int) comp_seq );
						const ScaleEntry volts = {"V" ,  0, 1.0};
						column.type = EngineConfig::TrajectoryLogger::LogColumn::Type::TOPLEVEL_STATE;
						column.value_type = EngineConfig::TrajectoryLogger::LogColumn::ValueType::F32;
						
						// TODO multinode - global state is not global anymore, redirect to comm buffers !
						size_t global_idx_V = tabs.global_state_f32_index[work_unit] + pig.GetVoltageStatevarIndex( comp_seq ); // LATER change for split cell?
						column.entry = global_idx_V;
						column.scaleFactor = Scales<Voltage>::native.ConvertTo(1, volts);
						// printf("\n\n\n record %zd \n\n\n", global_idx_V);
						break;
					}
					case Simulation::LemsQuantityPath::SegmentPath::Type::CALCIUM_INTRA:
					case Simulation::LemsQuantityPath::SegmentPath::Type::CALCIUM2_INTRA:
					{
						const ScaleEntry millimolar = {"mM" ,  0, 1.0}; // the reported scale is in fundamental units: mol/m^3 that is mmol/L
						column.type = EngineConfig::TrajectoryLogger::LogColumn::Type::TOPLEVEL_STATE;
						column.value_type = EngineConfig::TrajectoryLogger::LogColumn::ValueType::F32;
						
						const auto &comp_impl = pig.comp_implementations.at(comp_seq);
						const auto &comp_def  = pig.comp_definitions.at(comp_seq);
						
						Int Ca_seq = bioph.Ca_species_seq;
						const char *sCalcium = "calcium";
						if( path.segment.type == Simulation::LemsQuantityPath::SegmentPath::Type::CALCIUM2_INTRA ){
							Ca_seq = bioph.Ca2_species_seq;
							sCalcium = "calcium2";
						} 
						
						if( Ca_seq < 0 ){
							printf("internal error: logged biophysics missing %s", sCalcium);
							return false;
						}
						if( !comp_impl.concentration.count(Ca_seq) ){
							printf("internal error: logged biophysics missing %s impl", sCalcium);
							return false;
						}
						if( !comp_def.ions.count(Ca_seq) ){
							printf("internal error: logged biophysics missing %s def", sCalcium);
							return false;
						}
						
						const auto &calcimpl = comp_impl.concentration.at(Ca_seq);
						const auto &calcinst = comp_def.ions.at(Ca_seq);
						const auto &calcconc = conc_models.get(calcinst.conc_model_seq);
						
						ptrdiff_t Index_CaConcIn = calcimpl.Index_Intra;
						if( calcconc.type == ConcentrationModel::COMPONENT ){
							// get it from LEMS component
							
							Int comp_type_seq = calcconc.component.id_seq;
							if( comp_type_seq < 0 ){
								printf("internal error: lems quantity path for %s: missing component type\n", sCalcium);
								return false;
							}
							const ComponentType &comp_type = component_types.get(comp_type_seq);
							
							const auto exposure_seq = comp_type.common_exposures.concentration_intra;
							if( exposure_seq < 0 ){
								printf("internal error: lems quantity path for %s: missing component exposure %d\n", sCalcium, (int) exposure_seq);
								return false;
							}
							const auto &exposure = comp_type.exposures.get(exposure_seq);
							
							if( exposure.type == ComponentType::Exposure::STATE ){
								Index_CaConcIn = calcimpl.component.statevars_to_states[ exposure.seq ].index;
							}
							else{
								printf("error: lems quantity path for %s is not a state variable; this is not supported yet\n", sCalcium );
								return false; // perhaps fix LATER
							}
						}
						
						if( Index_CaConcIn < 0 ){
							printf("internal error: logged biophysics missing %s impl idx", sCalcium);
							return false;
						}
						
						// TODO multinode - global state is not global naymore, redirect to comm buffers !
						size_t global_idx_CaconcIn = tabs.global_state_f32_index[work_unit] + Index_CaConcIn; //TODO change for split cell?
						column.entry = global_idx_CaconcIn;
						column.scaleFactor = Scales<Concentration>::native.ConvertTo(1, millimolar);
						// printf("\n\n\n record %zd \n\n\n", global_idx_V);
						break;
					}
					default: {
						printf("column %ld for segment-located data writer %s not supported yet \n", col_seq, output_filepath.c_str());
						return false;
					}
				}
				
				break;
			}
			
			case Simulation::LemsQuantityPath::Type::CHANNEL: {
				
				if( !MustBePhysicalCell( cell_type ) ) return false;
				auto &pig = sig.physical_cell;
				
				switch(path.channel.type){
					
					case Simulation::LemsQuantityPath::ChannelPath::Type::Q: {
						
						column.type = EngineConfig::TrajectoryLogger::LogColumn::Type::TOPLEVEL_STATE;
						column.value_type = EngineConfig::TrajectoryLogger::LogColumn::ValueType::F32;
						
						// get states for the specific ion channel distribution
						const auto &comp_impl = pig.comp_implementations[comp_seq];
						size_t sig_Q_offset = comp_impl.channel[path.channel.distribution_seq].per_gate[path.channel.gate_seq].Index_Q;
						if(sig_Q_offset < 0){
							// TODO check for complex gates
							printf("column %ld for ion channel-located composite Q data writer %s not supported yet \n", col_seq, output_filepath.c_str());
							return false;
						}
						
						column.entry = tabs.global_state_f32_index[work_unit] + sig_Q_offset; //TODO change for split cell?
						column.scaleFactor = 1;
						// printf("\n\n\n record %zd \n\n\n", global_idx_V);
						break;
					}
					
					default: {
						printf("column %ld for ion channel-located data writer %s not supported yet \n", col_seq, output_filepath.c_str());
						return false;
					}
				}
				
				break;
			}
			
			case Simulation::LemsQuantityPath::Type::SYNAPSE:{
				printf("column %ld for data writer %s not supported yet : synapse path\n", col_seq, output_filepath.c_str());
				return false;
			}
			case Simulation::LemsQuantityPath::Type::INPUT:{
				printf("column %ld for data writer %s not supported yet : input path\n", col_seq, output_filepath.c_str());
				return false;
			}
			case Simulation::LemsQuantityPath::Type::CELL:{
				
				if( !MustBeArtificialCell( cell_type ) ) return false;
				auto &aig = sig.artificial_cell;
				
				Int comp_type_seq = cell_type.artificial.component.id_seq;
				if( comp_type_seq < 0 ){
					printf("internal error: lems quantity path for artificial cell: none native\n");
					return false;
				}
				
				const ComponentType &comp_type = component_types.get(comp_type_seq);
				const auto &namespace_thing_seq = path.cell.lems_quantity_path.namespace_thing_seq;
				const auto refer_thing = comp_type.name_space.get(namespace_thing_seq);
				
				ptrdiff_t Index_Statevar = -1;
				if( refer_thing.type == ComponentType::NamespaceThing::STATE ){
					Index_Statevar = aig.component.statevars_to_states[ refer_thing.seq ].index;
				}
				else{
					printf("error: lems quantity path for artificial cell is not a state variable; this is not supported yet\n");
					return false; // perhaps fix LATER
				}
				
				column.type = EngineConfig::TrajectoryLogger::LogColumn::Type::TOPLEVEL_STATE;
				column.value_type = EngineConfig::TrajectoryLogger::LogColumn::ValueType::F32;
				
				// TODO multinode - global state is not global anymore, redirect to comm buffers !
				size_t global_idx = tabs.global_state_f32_index[work_unit] + Index_Statevar;
				column.entry = global_idx;
				
				Dimension dim = comp_type.getNamespaceEntryDimension(namespace_thing_seq);
				// TODO scales should be strictly defined for lems quantities !					
				const LemsUnit native = dimensions.GetNative(dim);
				// const LemsUnit desired = native;					
				// use fundamental units for now
				const ScaleEntry fundamental_units = {"fundamental units" ,  0, 1.0}; 
				const LemsUnit desired = fundamental_units;
				
				column.scaleFactor = native.ConvertTo(1, desired);
				// printf("\n\n\n record %zd \n\n\n", global_idx_V);
				break;
			}
			
			case Simulation::LemsQuantityPath::Type::ION_POOL: // TODO
			default: {
				printf("column %ld for data writer %s not supported yet : cell-based path type %d\n", col_seq, output_filepath.c_str(), path.type);
				return false;
			}
			}
		}
		else{
			printf("column %ld for data writer %s not supported yet : non-cell-based path type %d \n", col_seq, output_filepath.c_str(), path.type);
			return false;
		}
		
		return true;
	};
	
	
	bool i_log_the_data = true;
	
	#ifdef USE_MPI
	if( my_mpi.rank == 0 ){
		
		i_log_the_data = true;
		// form the data structures, send requests for whatever is remote
		
	}
	else{
		i_log_the_data = false;
		// wait for remote requests from logger node to emerge
	}
	#endif
	
	if( i_log_the_data ){
		
		engine_config.trajectory_loggers.resize( sim.data_writers.contents.size() );
		for( Int daw_seq = 0; daw_seq < (Int)sim.data_writers.contents.size(); daw_seq++ ){
			
			auto &daw = sim.data_writers.get(daw_seq);
			
			EngineConfig::TrajectoryLogger &logger = engine_config.trajectory_loggers[daw_seq];
			
			logger.logfile_path = daw.fileName;
			
			for( Int col_seq = 0; col_seq < (Int)daw.output_columns.contents.size(); col_seq++ ){
				
				const auto &col = daw.output_columns.get(col_seq);
				const auto &path = col.quantity;
				
				EngineConfig::TrajectoryLogger::LogColumn column;
				
				if( !Implement_LoggerColumn( daw_seq, col_seq, daw.fileName, path, column ) ) return false;
				
				logger.columns.push_back(column);
			}
			
		}
	}
	//FIXME add event writers
	
	#ifdef USE_MPI
	
	printf("Determining recvlists...\n"); fflush(stdout);
	// Now, each node informs the others on what information streams it needs, in a symbolic(NeuroML-based) format
	if( config.debug_netcode ){
	Say("Recv");
	for( const auto &keyval : recv_lists ){
		Say("from node %d:", keyval.first);
		
		const auto &recv_list = keyval.second;
		for( const auto &keyval : recv_list.vpeer_refs ){
			const PointOnCellLocator &loc = keyval.first;
			std::string say_refs = "\tVpeer " + loc.toPresentableString() + " to remap refs: ";
			
			for( TabEntryRef_Packed ref : keyval.second ){
				say_refs += presentable_string(ref) + " ";
			}
			Say("%s", say_refs.c_str());
		}
		
		for( const auto &keyval : recv_list.spike_refs ){
			const PointOnCellLocator &loc = keyval.first;
			std::string say_refs = "\tSpikes " + loc.toPresentableString() + " to trigger refs: ";
			
			for( TabEntryRef_Packed ref : keyval.second ){
				say_refs += presentable_string(ref) + " ";
			}
			Say("%s", say_refs.c_str());
		}
		
		for( const auto &daw : recv_list.daw_refs ){
			std::string say_refs = "\tDaw " + daw.toPresentableString() + " to log ";
			
			Say("%s", say_refs.c_str());
		}
	}
	}
	
	printf("Exchanging recvlists...\n"); fflush(stdout);
	// I send recvlists to nodes, for them to send to me
	std::map< int, std::vector<char> > recvlists_encoded;
	
	// Nodes send recvlists to nodes, for me to send to them
	std::map< int, std::vector<char> > sendlists_encoded;
	
	// encode the recvlists for transmission
	for( const auto &keyval : recv_lists ){
		int other_rank = keyval.first;
		const RecvList &recvlist = keyval.second;
		std::string enc;
		// first, the header oine
		enc += accurate_string( recvlist.vpeer_refs.size() ) + " "
			+  accurate_string( recvlist.daw_refs  .size() ) + " "
			+  accurate_string( recvlist.spike_refs.size() ) + "\n" ;
		
		for( const auto &keyval : recvlist.vpeer_refs ){
			const auto &loc = keyval.first;
			loc.toEncodedString( enc );
			enc += "\n";
		}
		for( const DawRef &daw_ref : recvlist.daw_refs ){
			daw_ref.toEncodedString( enc );
			enc += "\n";
		}
		for( const auto &keyval : recvlist.spike_refs ){
			const auto &loc = keyval.first;
			loc.toEncodedString( enc );
			enc += "\n";
		}
		
		auto &recvlist_encoded = recvlists_encoded[other_rank];
		AppendToVector( recvlist_encoded, enc );
		recvlist_encoded.push_back('\0');
	}
	
	// dbg output encoded
	if( config.debug_netcode ){
	Say("Send Recvlist");
	for( const auto &keyval : recvlists_encoded ){
		int other_rank = keyval.first;
		const auto &recvlist_encoded = keyval.second;
		
		Say("to node %d, %s", other_rank, recvlist_encoded.data());
	}
	}
	// very much like Alltoallv, but communications are made with only existing connections (not the whole cartesian product)
	// this should really shine in large amounts of nodes (hundreds? or even tens?)
	// MPI_Graph_Neighbor* is similar, but it does not allow for automatic discovery of comm-adjacent nodes
	auto ExchangeLists = [  ]( const MPI_Datatype &datatype, const auto &sent_vectors_hashmap, auto & received_vectors_hashmap ){
		
		// skip, let them communicate
		
		const int TAG_LIST_SIZE = 1;
		const int TAG_LIST_LIST = 0;
		const int TAG_LIST_RECEIVED = 3;
		
		const bool use_Ssend = false; // for the msg list, so confiramtion comes through MPI itself
		
		// the steps are as following:
		// each node sends send lists to the nodes it receives from
		// and sets a "any" message sink to receive such lists from any receiver node it has to send to
		// the message sink is also used to receive confirmation of the send lists having been received !
		// 	this is important so that all send lists are sure to not be in the air, when listening stops
		// the process is completed by running a poll on whether all nodes have received confirmation that their send lists were accepted
		
		// an ordered sequence of the nodes I send recvlists to
		std::vector<int> send_seq_to_node;
		for( auto keyval : sent_vectors_hashmap ){
			send_seq_to_node.push_back(keyval.first);
		}
		int nSends = (int)send_seq_to_node.size();
		
		const int SEND_CODE_SUCCESS = 12345678;
		
		std::vector< MPI_Request > own_req_list;
		// Net receiver sends its recvlists to the nodes responsible to send the data
		
		// sends the sizes first (or more general headers, could be useful for indirect discovery too ?)
		int req_send_list_size_offset = own_req_list.size();
		own_req_list.insert( own_req_list.end(), nSends, MPI_REQUEST_NULL );
		
		// may send the lists too, or send them after connections are established (through list size headers)
		int req_send_list_offset = own_req_list.size();
		own_req_list.insert( own_req_list.end(), nSends, MPI_REQUEST_NULL );
		
		// and may get a confirmation from sender (or use Ssend to do this instead)
		int req_recv_list_confirm_offset = own_req_list.size();
		own_req_list.insert( own_req_list.end(), nSends, MPI_REQUEST_NULL );
		
		// These keep track of what the data receiver needs to send, and which sends are done
		std::vector<int> list_send_sizes(nSends,-1);	
		std::vector<int> recvlist_replies(nSends,-1);
		std::vector<bool> recvlist_replies_received( nSends, false );
		// keep a handy count of recvlist_replies_received
		int waiting_responses = 0;
		int waiting_responses_buffer = INT_MAX; // plus temp.storage for async polling
		
		
		// All nodes participate in the poll, to know when the connection discovery process is done
		int req_poll_offset = own_req_list.size();
		own_req_list.push_back( MPI_REQUEST_NULL );
		
		
		// Potential net sender needs to accept inbound messages from any other node
		int req_recv_list_size_offset = own_req_list.size();
		own_req_list.push_back( MPI_REQUEST_NULL );
		int recv_list_size_buffer = -1;
		
		// and has to track a dynamic, unknown number of inbound connections for the unknown-sized lists
		std::set<int> receiving_recvlist_from;
		
		std::map< int, MPI_Request > emergent_req_send_confirm;
		std::map< int, MPI_Request > emergent_req_recv_list; // for all messages, due to unknown size
		std::map< int, int > emergent_recv_list_size; // for all messages, due to unknown size
		
		
		// for Waitsome/Testsome
		std::vector< MPI_Status > status_buf;
		std::vector< int > status_buf_idx;
		
		// First, emit all messages
		
		// what if comm buffer is full? MPI runtime should then not put more messages, just keep them pending
		for( int i_recv = 0; i_recv < nSends; i_recv++ ){
			
			auto other_rank = send_seq_to_node[i_recv];
			const auto &sent_list = sent_vectors_hashmap.at(other_rank);
					
			if( sent_list.empty() ) continue;
			
			list_send_sizes[i_recv] = sent_list.size();
			
			// Sends to the same outbound node are alwyas ordered
			// send size msg first
			
			MPI_Isend( &list_send_sizes[i_recv], 1, MPI_INT, other_rank, TAG_LIST_SIZE, MPI_COMM_WORLD, &own_req_list[ req_send_list_size_offset + i_recv]);
			
			if(use_Ssend){
				// don't need to recv a confirmation message
				MPI_Issend( sent_list.data(), sent_list.size(), datatype, other_rank, TAG_LIST_LIST, MPI_COMM_WORLD, &own_req_list[ req_send_list_offset + i_recv]);
			}
			else{
				// recv a notifiation that the recvlist has been received 
				MPI_Isend( sent_list.data(), sent_list.size(), datatype, other_rank, TAG_LIST_LIST, MPI_COMM_WORLD, &own_req_list[ req_send_list_offset + i_recv]);
				
				MPI_Irecv( &recvlist_replies[i_recv], 1, MPI_INT, other_rank, TAG_LIST_RECEIVED, MPI_COMM_WORLD, &own_req_list[ req_recv_list_confirm_offset + i_recv]);
			}
			
			waiting_responses++;
		}
		
		auto Setup_Irecv_ListSize = [ & ]( ){
			MPI_Irecv( &recv_list_size_buffer, 1, MPI_INT, MPI_ANY_SOURCE, TAG_LIST_SIZE, MPI_COMM_WORLD, &own_req_list[ req_recv_list_size_offset ]);
		};
		Setup_Irecv_ListSize();
		
		// and get a random recv list too, for lists I should send to
		// MPI_Message probe_list_msg; // one for probe
		// int probe_ready;
		
		// Now, wait for incoming messages
		
		typedef std::chrono::high_resolution_clock::time_point Time;
		typedef std::chrono::duration<double> TimeDiff_Sec;
		Time tStart = std::chrono::high_resolution_clock::now();
		
		Time tLastPoll = tStart;
		auto getSecs = [  ]( Time tEnd, Time tstart ){
			return std::chrono::duration_cast< TimeDiff_Sec >(tEnd - tstart).count();
		};
		
		const double POLL_PERIOD_SEC = 0.1; // TODO make dynamic
		
		bool done = false;
		bool waiting_poll = false;
		int poll_result = -1;
		int poll_serial_no = 0;
		
		Say("Reqs %d %d %d %d %d %d", (int)own_req_list.size(), req_send_list_offset, req_send_list_size_offset, req_poll_offset, req_recv_list_confirm_offset, req_recv_list_size_offset);
		
		
		while( !done ){
			
			// don't log this when spinning
			// Say("Event loop start");
			
			Time tNow = std::chrono::high_resolution_clock::now();
			
			// NB: should spin anyway, because this is a time driven event -- not impossible for all communications to complete 
			if( getSecs( tNow, tLastPoll ) > POLL_PERIOD_SEC && !waiting_poll ){
				// start async poll
				Say("Start poll %d", poll_serial_no);
				waiting_responses_buffer = waiting_responses;
				MPI_Iallreduce( &waiting_responses_buffer, &poll_result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, &own_req_list[req_poll_offset] );
				poll_serial_no++;
				waiting_poll = true;
			}
			
			// Receive messages, participate in the polls, everything
			// Needs to spin because we also need to Recv an unknown-sized message -
			// better to spin than to waste another roundtrip, maybe?
			
			auto HandleRecvlistComplete = [ &recvlist_replies_received, &waiting_responses ]( int i_recv ){
				if( !recvlist_replies_received[i_recv] ){
					recvlist_replies_received[i_recv] = true;
					waiting_responses--;
				}
				else{
					assert(false); // how can a second mesage be sent ?
				}
			};
			
			
			// Poll for the static requests (recvlist-sender-side plus recvlist receiver for size/header)
			status_buf.resize( own_req_list.size() );
			status_buf_idx.resize( own_req_list.size() );
			int outcount;
			int ret = 0;
			
			ret = MPI_Testsome( own_req_list.size(), own_req_list.data(), &outcount, status_buf_idx.data(), status_buf.data() );
			
			bool success = ( ret == 0 );
			bool req_errors = ( ret == MPI_ERR_IN_STATUS );
			bool other_error = !( success || req_errors );
			
			if( req_errors ){
				Say("MPI POll Req Error, check it !");
				MPI_Abort( MPI_COMM_WORLD, 5 );
				abort();
			}
			if( other_error ){
				Say("MPI Poll Misc Error, check it !");
				MPI_Abort( MPI_COMM_WORLD, 5 );
				abort();				
			}
			
			for( int idx = 0; idx < outcount; idx++ ){
				int req_idx = status_buf_idx[idx];
				const MPI_Status &status = status_buf[idx];
				
				
				Say("Req %d finished", req_idx);
				
				if( req_send_list_offset <= req_idx && req_idx < req_send_list_offset + nSends ){
					int i_recv = req_idx - req_send_list_offset;
					
					Say("Send recvlist %d ( to node %d ) finished", i_recv, send_seq_to_node.at(i_recv));
					
					// send was done
					
					// if it's Ssend, send completion is as good as an expicit ack message !
					if( use_Ssend ){
						HandleRecvlistComplete(i_recv);
					}
				}
				else if( req_send_list_size_offset <= req_idx && req_idx < req_send_list_size_offset + nSends ){
					int i_recv = req_idx - req_send_list_size_offset;
					
					Say("Send recvlist size %d ( to node %d ) finished", i_recv, send_seq_to_node.at(i_recv));
					
					// nothing meaningful to do
				}
				else if( req_idx == req_poll_offset ){
					Say("Received poll result of %d", poll_result);
					tLastPoll = tNow;
					waiting_poll = false;
					
					if( poll_result == 0 ){
						// finished, yay!!
						done = true;
						// maybe don't leave yet, just in case there is cleanup to be done in this loop
					}
					
					// otherwise, wait till next poll
				}
				else if( req_recv_list_confirm_offset <= req_idx && req_idx < req_recv_list_confirm_offset + nSends ){
					int i_recv = req_idx - req_recv_list_confirm_offset;
					
					int other_rank = status.MPI_SOURCE;
					auto should_be_other_rank = send_seq_to_node[i_recv];
					
					int recv_msg = recvlist_replies[i_recv] ;
					
					Say("Received confirmation from node %d, is %d", other_rank, recv_msg );
					
					assert( recv_msg == SEND_CODE_SUCCESS );
					assert( other_rank == should_be_other_rank ); (void)should_be_other_rank;
					
					HandleRecvlistComplete(i_recv);
				}
				else if( req_idx == req_recv_list_size_offset ){
					
					int other_rank = status.MPI_SOURCE;
					int tag = status.MPI_TAG;
					
					int sendlist_size = recv_list_size_buffer;
					
					Say("Received recvlist size from %d, length %d, yag %d", other_rank, sendlist_size, tag );
					
					if( tag != TAG_LIST_SIZE ){
						assert(false);
					}
					
					if( receiving_recvlist_from.count(other_rank) ){
						Say("But already received from that node !!");
						continue;
						assert(false);
					}
					
					
					receiving_recvlist_from.insert( other_rank );
					emergent_recv_list_size[other_rank] = sendlist_size;
					
					// receive the lists too, in this phase
					// stay in the poll by not sending the confirmation yet
					
					auto Recv_Recvlist = [ &received_vectors_hashmap, &emergent_req_recv_list, &datatype ]( int other_rank, int sendlist_size ){
						// allocate a send list, directly
						// to get the received_list buffer
						auto &received_list = received_vectors_hashmap[other_rank];
						
						received_list.resize(sendlist_size);
						
						// also emergent entries
						auto &req = emergent_req_recv_list[other_rank];
						
						MPI_Irecv( received_list.data(), received_list.size(), datatype, other_rank, TAG_LIST_LIST, MPI_COMM_WORLD, &req );
						
						Say("Receiving recvlist from %d", other_rank);
					}; 
					Recv_Recvlist( other_rank, sendlist_size );
					
					
					// Also, since the recv completed, set up the recv for the next prospective connection
					Setup_Irecv_ListSize();
				}
				else{
					// why??
					assert(false);
				}
			}
			
			// But also check from the recvlist-receiver-side (schedule)
			for( auto &keyval : emergent_req_recv_list ){
				const auto &other_rank = keyval.first;
				auto &req_recv_list = keyval.second;
				
				if( req_recv_list == MPI_REQUEST_NULL ) continue; // TODO eras finished ones, somehow
				
				int flag = 0;
				MPI_Status status;
			
				MPI_Test( &req_recv_list, &flag, &status);
				if( flag ){
					Say("Received recvlist from %d", other_rank);
					// reception complete
					
					// and send a confirmation
					auto Send_Recvlist_Confirmation = [ &SEND_CODE_SUCCESS, &use_Ssend, &emergent_req_send_confirm ]( int other_rank ){
						if( use_Ssend ){
							// sender is implicitly notified
						}
						else{
							MPI_Isend( &SEND_CODE_SUCCESS, 1, MPI_INT, other_rank, TAG_LIST_RECEIVED, MPI_COMM_WORLD, &emergent_req_send_confirm[ other_rank ]);
						};
					};
					Send_Recvlist_Confirmation( other_rank );
					
					// could do something inline with received lists, perhaps a functor LATER
				}
				else{
					// not ready, pass
				}
			}
			
			
			// perhaps yield some time CPU time/energy, in a cooperative application LATER
		
		}
		
		// cancel the persistent random-origin recv, and cleanup all requests
		MPI_Cancel( &own_req_list[req_recv_list_size_offset] );
		
		// the rest is probably not necessary because no message should be flying, still cleanup
		MPI_Waitall( own_req_list.size(), own_req_list.data(), MPI_STATUSES_IGNORE);
		// and emergent sends, too
		for( auto &keyval : emergent_req_send_confirm ){
			auto &req = keyval.second;
			MPI_Wait( &req, MPI_STATUS_IGNORE );
		}
		for( auto &keyval : emergent_req_recv_list ){
			auto &req = keyval.second;
			MPI_Wait( &req, MPI_STATUS_IGNORE );
		}
		
		
		Time tEnd = std::chrono::high_resolution_clock::now();
		Say("Finished exchanging send lists in %f sec.", getSecs( tEnd, tStart ) );
		
		// assert(false);
		
		// and form the mirrors
	};
	
	ExchangeLists( MPI_CHAR, recvlists_encoded, sendlists_encoded );
	
	// dbg output encoded
	if( config.debug_netcode ){
	Say("Received Recvlists");
	for( const auto &keyval : sendlists_encoded ){
		int other_rank = keyval.first;
		const auto &recvlist_encoded = keyval.second;
		
		Say("from node %d, %s", other_rank, recvlist_encoded.data());
	}
	}
	
	// decode the recvlists
	for( auto &keyval : sendlists_encoded ){
		
		auto other_rank = keyval.first;
		auto &enc = keyval.second; // newlines will be replaced with nulls to ease parsing
		
		auto &sendlist = send_lists[other_rank];
		assert( enc.size() );
		// convert enc to null-terminated lines, for convenience
		std::vector<char *> lines;
		lines.push_back( enc.data() );
		for( size_t i = 0; i < enc.size() ; i++ ){
			if( enc[i] == '\n' ){
				enc[i] = '\0';
				if( i + 1 < enc.size() ) lines.push_back( enc.data() + i + 1 );
			}
		}
		
		if( config.debug_netcode ){
		Say("Lines: %zd", lines.size() );
		for( int i = 0; i < (int)lines.size() ; i++ ){
			Say("%d:\t%s", i, lines[i]);
			Say("end-----\n\n");
		}
		}
		
		Int vpeers, daws, spikes;
		sscanf(lines[0], "%ld %ld %ld", &vpeers, &daws, &spikes );
		if( config.debug_netcode ){
			Say("%ld %ld %ld <- %s", vpeers, daws, spikes, lines[0]);
		}
		Int vpeer_idx = 1;
		Int daw_idx = vpeer_idx + vpeers;
		Int spike_idx = daw_idx + daws;
		
		
		sendlist.vpeer_sources.resize(vpeers);
		for(int i = 0; i < vpeers; i++){
			auto ret = sendlist.vpeer_sources[i].fromEncodedString( lines[ vpeer_idx + i ] );
			if(!ret) Say( "fail %s", lines[ vpeer_idx + i ] );
			assert(ret);
		}
		sendlist.daw_refs.resize(daws);
		for(int i = 0; i < daws; i++){
			auto ret = sendlist.daw_refs[i].fromEncodedString( lines[ daw_idx + i ] );
			if(!ret) Say( "fail %s", lines[ daw_idx + i ] );
			assert(ret);
		}
		sendlist.spike_sources.resize(spikes);
		for(int i = 0; i < spikes; i++){
			auto ret = sendlist.spike_sources[i].fromEncodedString( lines[ spike_idx + i ] );
			if(!ret) Say( "fail %s", lines[ spike_idx + i ] );
			assert(ret);
		}
	}
	
	// dbg output encoded symbolic
	if( config.debug_netcode ){
	Say("Send");
	for( const auto &keyval : send_lists ){
		Say("to node %d:", keyval.first);
		
		const auto &send_list = keyval.second;
		for( const auto &loc : send_list.vpeer_sources ){
			std::string say_refs = "\tVpeer " + loc.toPresentableString() ;
			Say("%s", say_refs.c_str());
		}
		
		for( const auto &loc : send_list.spike_sources ){
			std::string say_refs = "\tSpikes " + loc.toPresentableString() ;
			Say("%s", say_refs.c_str());
		}
		
		for( const auto &daw : send_list.daw_refs ){
			std::string say_refs = "\tDaw " + daw.toPresentableString() ;
			Say("%s", say_refs.c_str());
		}
	}
	}
	
	// now construct the send and recv mirrors, and remap
	
	// construct and remap for send_lists
	for( const auto &keyval : send_lists ){
		auto other_rank = keyval.first;
		const auto &send_list = keyval.second;
		
		auto &send_list_impl = engine_config.sendlist_impls[other_rank];
		
		send_list_impl.vpeer_positions_in_globstate.resize( send_list.vpeer_sources.size() );
		for( size_t i = 0; i < send_list.vpeer_sources.size() ; i++ ){
			const auto &loc = send_list.vpeer_sources.at(i);
			send_list_impl.vpeer_positions_in_globstate[i] = GetCompartmentVoltageStatevarIndex_Global( loc );
		}
		
		send_list_impl.daw_columns.resize( send_list.daw_refs.size() );
		// also realize any data logger columns that access values local to this node
		for( size_t i = 0; i < send_list.daw_refs.size() ; i++ ){
			const auto &ref = send_list.daw_refs.at(i);
			
			// NB make sure this node knows about this daw
			const auto &daw = sim.data_writers.get(ref.daw_seq);
			const auto &col = daw.output_columns.get(ref.col_seq);
			const auto &path = col.quantity;
			
			auto &impl_col = send_list_impl.daw_columns[i];
			
			if( !Implement_LoggerColumn( ref.daw_seq, ref.col_seq, daw.fileName, path, impl_col) ) return false;
		}
		
		// allocate mirror buffers for spike triggers
		send_list_impl.spike_mirror_buffer = tabs.global_tables_state_i64_arrays.size();
		tabs.                                     global_tables_state_i64_arrays.emplace_back();
		auto &tab = tabs.                         global_tables_state_i64_arrays.back();
		// TODO pack the boolean vectors
		tab.resize( send_list.spike_sources.size(), 0);
		// and add extra notification entries to the spike sources
		for( size_t i = 0; i < send_list.spike_sources.size() ; i++ ){
			const auto &loc = send_list.spike_sources.at(i);
			
			size_t global_idx_T_spiker;
			if( !GetCompartmentSpikerImplementation_Global( loc, global_idx_T_spiker) ) return false;
			auto packed_id = GetEncodedTableEntryId( send_list_impl.spike_mirror_buffer, i );
			
			tabs.global_tables_const_i64_arrays[global_idx_T_spiker].push_back( packed_id );
		}
	}
	// The final send buffer for these will be allocated at run time
	
	// construct and remap for recv_lists
	for( const auto &keyval : recv_lists ){
		auto other_rank = keyval.first;
		const auto &recv_list = keyval.second;
		
		auto &recv_list_impl = engine_config.recvlist_impls[other_rank];
		
		// allocate mirror buffers for values continuously being sent
		recv_list_impl.value_mirror_size = recv_list.vpeer_refs.size() + recv_list.daw_refs.size();
		recv_list_impl.value_mirror_buffer = tabs.global_tables_state_f32_arrays.size();
		tabs.                                     global_tables_state_f32_arrays.emplace_back();
		auto &value_mirror = tabs.                global_tables_state_f32_arrays.back();
		value_mirror.resize( recv_list_impl.value_mirror_size, 5555); // XXX set as NAN after debug
		for( int i = recv_list.vpeer_refs.size(); i < (int)value_mirror.size(); i++ ) value_mirror[i] = 4444;
		// NB walk through these requested values, in the same order they were requested in the recv list
		long long value_mirror_table = recv_list_impl.value_mirror_buffer;
		int value_mirror_entry = 0;
		for( const auto &keyval : recv_list.vpeer_refs ){
			const auto & remap_ref_list = keyval.second;
			
			for( TabEntryRef_Packed ref_packed : remap_ref_list ){
				TabEntryRef ref = GetDecodedTableEntryId(ref_packed);
				TabEntryRef_Packed remapped_ref = GetEncodedTableEntryId( value_mirror_table, value_mirror_entry );
				// printf("reff %llx -> %lld %d -> %llx, %zx %zx\n", ref_packed, ref.table, ref.entry, remapped_ref, tabs.global_tables_const_i64_arrays.size() , tabs.global_tables_const_i64_arrays[ref.table].size());
				tabs.global_tables_const_i64_arrays[ref.table][ref.entry] = remapped_ref;
			}
			
			value_mirror_entry++;
		}
		
		
		// right next, the daw values
		for( const auto &daw_ref : recv_list.daw_refs ){
			
			assert(my_mpi.rank == 0);
			
			engine_config.trajectory_loggers[daw_ref.daw_seq].columns[daw_ref.col_seq].entry = value_mirror_entry;
			
			value_mirror_entry++;
		}
		
		// and also track where the spikes should be sent to
		int spike_mirror_entry = 0;
		recv_list_impl.spike_destinations.resize( recv_list.spike_refs.size() );
		for( const auto &keyval : recv_list.spike_refs ){
			const auto &ref_list = keyval.second;
			
			recv_list_impl.spike_destinations[spike_mirror_entry] = ref_list;
			
			spike_mirror_entry++;
		}
		#ifdef LOLOL
		#endif
	}
	
	// MPI_Finalize();
	// exit(1);
	#endif
	// some final info
	
	engine_config.work_items = tabs.callbacks.size(); // kind of obvious in hindsight
	engine_config.t_initial = 0; // might change LATER
	engine_config.t_final = engine_config.t_initial + sim.length;
	engine_config.dt = sim.step;
	
	// yay!
	printf("instantiation complete!\n");
	
	return true;
}

int main(int argc, char **argv){
	
	#ifdef USE_MPI
	// first of first of all, replace argc and argv
	// Modern implementations may keep MPI args from appearing anyway; non-modern ones still need this
	MPI_Init(&argc, &argv);
	#endif
	
	// first of all, set stdout,stderr to Unbuffered, for live output
	// this action must happen before any output is written !
	setvbuf(stdout, NULL, _IONBF, 0);
	setvbuf(stderr, NULL, _IONBF, 0);
	
	//sim parameters
	SimulatorConfig config;
	
	//print logo ^_^
	printf("--- Extensible Dynamics Engine for Networks ---\n");
	#ifndef BUILD_STAMP
	#define BUILD_STAMP __DATE__
	#endif
	printf("Build version " BUILD_STAMP "\n");
	
	#ifdef USE_MPI
	// Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &my_mpi.world_size);
    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi.rank);
	
	char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
	
	
	if( 1 || config.verbose){
		printf("Hello from processor %s, rank %d out of %d processors\n", processor_name, my_mpi.rank, my_mpi.world_size);
		
		if( my_mpi.rank != 0 ){
			char tmps[555];
			sprintf(tmps, "log_node_%d.gen.txt", my_mpi.rank);
			freopen(tmps,"w",stdout);
			stderr = stdout;
		}
	}
	
	#endif
	
	RunMetaData metadata;
	
	Model model; // TODO move to SimulatorConfig
	bool model_selected = false;
	
	timeval config_start, config_end;
	gettimeofday(&config_start, NULL);
	//-------------------> Read config options for run
	for (int i = 1; i < argc; i++){
		std::string arg(argv[i]);
		
		// model options
		if(arg == "nml"){
			//read JSON
			if(i == argc - 1){
				printf("cmdline: NeuroML filename missing\n");
				exit(1);
			}
			
			timeval nml_start, nml_end;
			gettimeofday(&nml_start, NULL);
			if(!( ReadNeuroML( argv[i+1], model, true, config.verbose) )){
				printf("cmdline: could not make sense of NeuroML file\n");
				exit(1);
			}
			gettimeofday(&nml_end, NULL);
			printf("cmdline: Parsed %s in %lf seconds\n", argv[i+1], TimevalDeltaSec(nml_start, nml_end));
			
			model_selected = true;
			
			i++;
		}
		// model overrides
		else if(arg == "rng_seed"){
			if(i == argc - 1){
				printf( "cmdline: %s value missing\n", arg.c_str() );
				exit(1);
			}
			const std::string sSeed = argv[i+1];
			long long seed;
			if( sscanf( sSeed.c_str(), "%lld", &seed ) == 1 ){
				config.override_random_seed = true;
				config.override_random_seed_value = seed;
			}
			else{
				printf( "cmdline: %s must be a reasonably-sized integer, not %s\n", arg.c_str(), sSeed.c_str() );
				exit(1);
			}
			
			i++; // used following token too
		}
		// solver options
		else if(arg == "cable_solver"){
			//read JSON
			if(i == argc - 1){
				printf( "cmdline: %s type missing\n", arg.c_str() );
				exit(1);
			}
			const std::string soltype = argv[i+1];
			if( soltype == "fwd_euler" ){
				printf("Cable solver set to Forward Euler: make sure system is stable\n");
				config.cable_solver = SimulatorConfig::CABLE_FWD_EULER;
			}
			else if( soltype == "bwd_euler" ){
				config.cable_solver = SimulatorConfig::CABLE_BWD_EULER;
			}
			else if( soltype == "auto" ){
				config.cable_solver = SimulatorConfig::CABLE_SOLVER_AUTO;
			}
			else{
				printf( "cmdline: unknown %s type %s; "
					"choices are auto, fwd_euler, bwd_euler\n", arg.c_str(),  soltype.c_str() );
				exit(1);
			}
			
			i++; // used following token too
		}
		// debugging options
		else if(arg == "verbose"){
			config.verbose = true;
		}
		else if(arg == "full_dump"){
			config.dump_raw_state_scalar = true;
			config.dump_raw_state_table = true;
			config.dump_raw_layout = true;
		}
		else if(arg == "dump_state_scalar"){
			config.dump_raw_state_scalar = true;
		}
		else if(arg == "dump_raw_layout"){
			config.dump_raw_layout = true;
		}
		else if(arg == "debug"){
			config.debug = true;
			config.debug_netcode = true;
		}
		else if(arg == "debug_netcode"){
			config.debug_netcode = true;
		}
		else if(arg == "-S"){
			config.output_assembly = true;
		}
		// optimization, compiler options
		else if(arg == "icc"){
			config.use_icc = true;
		}
		else if(arg == "gcc"){
			config.use_icc = false;
		}
		
		else{
			//unknown, skip it
			printf("cmdline: skipping unknown token \"%s\"\n", argv[i]);
		}
	}
	gettimeofday(&config_end, NULL);
	
	// handle model missing
	if(!model_selected){
		printf("error: NeuroML model not selected (select one with nml <file> in command line)\n");
		return 2;
	}
	
	metadata.config_time_sec = TimevalDeltaSec(config_start, config_end);
	
	//-------------------> create data structures and code based on the model
	timeval init_start, init_end;
	gettimeofday(&init_start, NULL);
	
	printf("Initializing model...\n");
	EngineConfig engine_config;
	RawTables tabs;
	if(!GenerateModel(model, config, engine_config, tabs)){
		printf("NeuroML model could not be created\n");
		exit(1);
	}
	
	//-------------------> crunch the numbers
	// set up printing in logfiles
	const int column_width = 16;
	char tmps_column[ column_width + 5 ];
	FixedWidthNumberPrinter column_fmt(column_width, '\t', 0);
	
	// open the logs, one for each logger
	std::vector<FILE *> trajectory_open_files;
	for(auto logger : engine_config.trajectory_loggers){
		#ifdef USE_MPI
		assert( my_mpi.rank == 0);
		#endif
		const char *path = logger.logfile_path.c_str();
		FILE *fout = fopen( path, "wt");
		if(!fout){
			auto errcode = errno;// NB: keep errno right away before it's overwritten
			printf("Could not open trajectory log \"%s\" : %s\n", path, strerror(errcode) );
			exit(1);
		}
		trajectory_open_files.push_back(fout);
	}
	
	// prepare engine for crunching
	
	printf("Allocating state buffers...\n");
	// allocate at least two state vectors, to iterate in parallel
	RawTables::Table_F32 state_one = tabs.global_initial_state; // eliminate redundancy LATER
	RawTables::Table_F32 state_two(tabs.global_initial_state.size(), NAN);
	
	auto	 		tables_state_f32_one = tabs.global_tables_state_f32_arrays;
	decltype(		tables_state_f32_one )  			tables_state_f32_two;
					tables_state_f32_two.reserve(		tables_state_f32_one.size());
	for( auto tab : tables_state_f32_one )      		tables_state_f32_two.emplace_back( tab.size(), NAN );
	
	auto     		tables_state_i64_one = tabs.global_tables_state_i64_arrays;
	decltype(		tables_state_i64_one)              	tables_state_i64_two;
					tables_state_i64_two.reserve(      	tables_state_i64_one.size() );
	for( auto tab : tables_state_i64_one )      		tables_state_i64_two.emplace_back( tab.size(), 0 );
	// now things need to be done a little differently, since for example trigger(and lazy?) variables of Next ought to be zero for results to make sense
	
	float *global_state_now = &state_one[0];
	float *global_state_next = &state_two[0];
	
	//also allocate pointer and size vectors, to use instead of silly std::vectors
	auto GetSizePtrTables = []( auto &tablist, auto &pointers, auto &sizes ){
		
		pointers.resize( tablist.size() );
		sizes.resize( tablist.size() );
		for(size_t i = 0; i < tablist.size(); i++){
			pointers[i] = tablist.at(i).data();
			sizes[i] = (long long)tablist.at(i).size();
		}
	};
	std::vector <long long> global_tables_const_f32_sizes; std::vector <Table_F32> global_tables_const_f32_arrays;
	std::vector <long long> global_tables_const_i64_sizes; std::vector <Table_I64> global_tables_const_i64_arrays;
	
	std::vector <long long> global_tables_state_f32_sizes; std::vector <Table_F32> global_tables_stateOne_f32_arrays; std::vector <Table_F32> global_tables_stateTwo_f32_arrays;
	std::vector <long long> global_tables_state_i64_sizes; std::vector <Table_I64> global_tables_stateOne_i64_arrays; std::vector <Table_I64> global_tables_stateTwo_i64_arrays;
	
	GetSizePtrTables(tabs.global_tables_const_f32_arrays, global_tables_const_f32_arrays, global_tables_const_f32_sizes);
	GetSizePtrTables(tabs.global_tables_const_i64_arrays, global_tables_const_i64_arrays, global_tables_const_i64_sizes);
	
	GetSizePtrTables(tables_state_f32_one, global_tables_stateOne_f32_arrays, global_tables_state_f32_sizes);
	GetSizePtrTables(tables_state_i64_one, global_tables_stateOne_i64_arrays, global_tables_state_i64_sizes);
	GetSizePtrTables(tables_state_f32_two, global_tables_stateTwo_f32_arrays, global_tables_state_f32_sizes);
	GetSizePtrTables(tables_state_i64_two, global_tables_stateTwo_i64_arrays, global_tables_state_i64_sizes);
	
	// also, set up the references to the flat vectors
	global_tables_const_f32_arrays[tabs.global_const_tabref] = tabs.global_constants.data();
	global_tables_const_f32_sizes [tabs.global_const_tabref] = tabs.global_constants.size();
	
	global_tables_stateOne_f32_arrays[tabs.global_state_tabref] = state_one.data();
	global_tables_stateTwo_f32_arrays[tabs.global_state_tabref] = state_two.data();
	global_tables_state_f32_sizes    [tabs.global_state_tabref] = state_one.size();
	
	Table_F32 *global_tables_stateNow_f32  = global_tables_stateOne_f32_arrays.data();
	Table_I64 *global_tables_stateNow_i64  = global_tables_stateOne_i64_arrays.data();
	Table_F32 *global_tables_stateNext_f32 = global_tables_stateTwo_f32_arrays.data();
	Table_I64 *global_tables_stateNext_i64 = global_tables_stateTwo_i64_arrays.data();
	
	gettimeofday(&init_end, NULL);
	metadata.init_time_sec = TimevalDeltaSec(init_start, init_end);
	
	if(config.dump_raw_layout){
		
		printf("Constants:\n");
		for(auto val : tabs.global_constants ) printf("%g \t", val); 
		printf("\n");
		
		printf("ConstIdx:\n");
		for(auto val : tabs.global_const_f32_index ) printf("%lld \t", val);
		printf("\n");
		printf("StateIdx:\n");
		for(auto val : tabs.global_state_f32_index ) printf("%lld \t", val);
		printf("\n");
		
		auto PrintTables = [](const auto &index, const auto &arrays){
			size_t next_tabchunk = 0;
			for(size_t i = 0; i < arrays.size() ;i++){
				
				if(next_tabchunk < index.size() && i == (size_t)index.at(next_tabchunk) ){
					printf("%zd", i);
					while( next_tabchunk < index.size() && i == (size_t)index.at(next_tabchunk) ) next_tabchunk++;
				}
				printf(" \t");
				printf(" %16p \t", arrays.at(i).data());
				for(auto val : arrays.at(i) ) printf("%s \t", (presentable_string(val)).c_str());
				printf("\n");
			}
		} ;
		auto PrintRawTables = [](const auto &index, const auto &arrays, const auto &sizes){
			size_t next_tabchunk = 0;
			for(size_t i = 0; i < arrays.size() ;i++){
				
				if(next_tabchunk < index.size() && i == (size_t)index.at(next_tabchunk) ){
					printf("%zd", i);
					while( next_tabchunk < index.size() && i == (size_t)index.at(next_tabchunk) ) next_tabchunk++;
				}
				printf(" \t");
				printf(" %16p \t", arrays.at(i));
				for( std::size_t j = 0; j < (size_t)sizes.at(i); j++ ){
					auto val = arrays.at(i)[j];
					printf("%s \t", (presentable_string(val)).c_str());
				}
				printf("\n");
			}
		} ;
		
		
		printf("TabConstF32: %zd %zd\n", tabs.global_table_const_f32_index.size(), tabs.global_tables_const_f32_arrays.size() );
		PrintTables(tabs.global_table_const_f32_index, tabs.global_tables_const_f32_arrays);
		printf("TabConstI64: %zd %zd\n", tabs.global_table_const_i64_index.size(), tabs.global_tables_const_i64_arrays.size() );
		PrintTables(tabs.global_table_const_i64_index, tabs.global_tables_const_i64_arrays);
		printf("TabStateF32: %zd %zd\n", tabs.global_table_state_f32_index.size(), tabs.global_tables_state_f32_arrays.size() );
		PrintTables(tabs.global_table_state_f32_index, tabs.global_tables_state_f32_arrays);
		printf("TabStateI64: %zd %zd\n", tabs.global_table_state_i64_index.size(), tabs.global_tables_state_i64_arrays.size() );
		PrintTables(tabs.global_table_state_i64_index, tabs.global_tables_state_i64_arrays);
		
		printf("RawStateI64:\n");
		PrintRawTables(tabs.global_table_state_i64_index, global_tables_stateOne_i64_arrays, global_tables_state_i64_sizes);
		
		printf("CallIdx:\n");
		for(auto val : tabs.callbacks ) printf("%p \t", val);
		printf("\n");
		
		
		printf("Initial state:\n");
		printf("TabStateOneF32:\n");
		PrintTables(tabs.global_table_state_f32_index, tables_state_f32_one);
		printf("TabStateOneI64:\n");
		PrintTables(tabs.global_table_state_i64_index, tables_state_i64_one);
		printf("TabStateTwoF32:\n");
		PrintTables(tabs.global_table_state_f32_index, tables_state_f32_two);
		printf("TabStateTwoI64:\n");
		PrintTables(tabs.global_table_state_i64_index, tables_state_i64_two);
		printf("Initial scalar state:\n");
		for(auto val : tabs.global_initial_state ) printf("%g \t", val);
		printf("\n");
		
	}
	// MPI_Finalize();
	// exit(1);
	
	#ifdef USE_MPI
	
	printf("Allocating comm buffers...\n");
	typedef std::vector<float> SendRecvBuf;
	std::vector<int> send_off_to_node;
	std::vector< SendRecvBuf > send_bufs;
	std::vector<int> recv_off_to_node;
	std::vector< SendRecvBuf > recv_bufs;
	
	for( const auto &keyval : engine_config.sendlist_impls ){
		send_off_to_node.push_back( keyval.first );
		send_bufs.emplace_back();
		// allocate as they come, why not
	}
	
	for( const auto &keyval : engine_config.recvlist_impls ){
		recv_off_to_node.push_back( keyval.first );
		recv_bufs.emplace_back();
		// allocate as they come, why not
	}
	std::vector<MPI_Request> send_requests( send_off_to_node.size(), MPI_REQUEST_NULL );
	std::vector<MPI_Request> recv_requests( recv_off_to_node.size(), MPI_REQUEST_NULL );
	
	// recv's have to be probed before recv'ing
	std::vector<bool> received_probes( recv_off_to_node.size(), false);	
	std::vector<bool> received_sends( recv_off_to_node.size(), false);
	
	#endif
	
	printf("Starting simulation loop...\n");
	
	// NB this is a simple tweak to avoid spawning more threads than work items.
	// TODO allow for diagnostic override, and amore flexibility since multithreading can slao be used for logging etc.
	// (how much time does it take to change the number of threads?)
	// NB: omp_get_max_threads actually refers to the very real nthreads-var that is also picked through the OMP_NUM_THREADS envvar.
	// That avriable is actually set through omp_set_num_threads
	if( engine_config.work_items < omp_get_max_threads() ){
		omp_set_num_threads( (int)engine_config.work_items );
	}
	
	// perform the crunching
	timeval run_start, run_end;
	gettimeofday(&run_start, NULL);
	
	double time = engine_config.t_initial;
	// need multiple initialization steps, to make sure the dependency chains of all state variables are resolved
	for( long long step = -3; time <= engine_config.t_final; step++ ){
		
		bool initializing = step <= 0;
		
		#ifdef USE_MPI
		auto NetMessage_ToString = []( size_t buf_value_len, const auto &buf ){
			std::string str;
			for( size_t i = 0; i < buf_value_len; i++ ){
				str += presentable_string( buf[i] ) + " ";
			}
			str += "| ";
			for( size_t i = buf_value_len; i < buf.size(); i++ ){
				str += presentable_string( EncodeF32ToI32( buf[i] ) ) + " ";
			}
			return str;
		};
		// Send info needed by other nodes
		// TODO try parallelizing buffer fill, see if it improves latency
		for( size_t idx = 0; idx < send_off_to_node.size(); idx++ ){
			auto other_rank = send_off_to_node.at(idx);
			const auto &sendlist_impl = engine_config.sendlist_impls.at(other_rank);
			
			auto &buf = send_bufs[idx];
			auto &req = send_requests[idx];
			
			// get the continuous_time values
			size_t vpeer_buf_idx = 0;
			size_t vpeer_buf_len = sendlist_impl.vpeer_positions_in_globstate.size();
			
			size_t daw_buf_idx = vpeer_buf_idx + vpeer_buf_len;
			size_t daw_buf_len = sendlist_impl.daw_columns.size();
			
			size_t buf_value_len = vpeer_buf_len + daw_buf_len;
			
			buf.resize(buf_value_len); // std::vector won't reallocate due to size reduction
			// otherwise implement appending the spikes manually
			
			// NB make sure these buffers are synchronized with CPU memory LATER
			for( size_t i = 0; i < sendlist_impl.vpeer_positions_in_globstate.size(); i++ ){
				size_t off = sendlist_impl.vpeer_positions_in_globstate[i];
				
				buf[ vpeer_buf_idx + i ] = global_state_now[off];
			}
			
			for( size_t i = 0; i < sendlist_impl.daw_columns.size(); i++ ){
				assert( my_mpi.rank != 0 && other_rank == 0 );
				auto &col = sendlist_impl.daw_columns[i];
				size_t off = col.entry;
				// also apply scaling, so receiving node won't bother
				
				buf[ daw_buf_idx + i ] = global_state_now[ off ] * col.scaleFactor ;
			}
			
			size_t spikebuf_off = sendlist_impl.spike_mirror_buffer;
			// get the spikes into the buffer (variable size)
			Table_I64 SpikeTable = global_tables_stateNow_i64[spikebuf_off];
			long long SpikeTable_size = global_tables_state_i64_sizes[spikebuf_off];
			
			for( int i = 0; i < SpikeTable_size; i++ ){
			
				// TODO packed bool buffers
				if( SpikeTable[i] ){
					// add index	
					buf.push_back(  EncodeI32ToF32(i) );
					// clear trigger flag for the timestep after the next one
					SpikeTable[i] = 0;
				}
			}
			if( config.debug_netcode ){
				Say("Send %d : %s", other_rank, NetMessage_ToString( buf_value_len, buf).c_str());
			}
			MPI_Isend( buf.data(), buf.size(), MPI_FLOAT, other_rank, MYMPI_TAG_BUF_SEND, MPI_COMM_WORLD, &req );
		}
		
		// Recv info needed by this node
		auto PostRecv = []( int other_rank, std::vector<float> &buf, MPI_Request &recv_req ){
			MPI_Irecv( buf.data(), buf.size(), MPI_FLOAT, other_rank, MYMPI_TAG_BUF_SEND, MPI_COMM_WORLD, &recv_req );
		};
		auto ReceiveList = [ &engine_config, &global_tables_stateNow_f32, &global_tables_stateNow_i64 ]( const EngineConfig::RecvList_Impl &recvlist_impl, std::vector<float> &buf ){
			
			// copy the continuous-time values
			size_t value_buf_idx = 0;
			float *value_buf = global_tables_stateNow_f32[ recvlist_impl.value_mirror_buffer ];
			// NB make sure these buffers are synchronized with CPU memory LATER
			for( ptrdiff_t i = 0; i < recvlist_impl.value_mirror_size; i++ ){
				value_buf[i] = buf[ value_buf_idx + i ];
				
				// global_state_now[ off + i] = buf[ value_buf_idx + i ];
			}
			
			// and deliver the spikes to trigger buffers
			for( int i = recvlist_impl.value_mirror_size; i < (int)buf.size(); i++ ){
				int spike_pos = EncodeF32ToI32( buf[i] );
				for( auto tabent_packed : recvlist_impl.spike_destinations[spike_pos] ){
					auto tabent = GetDecodedTableEntryId( tabent_packed );
					// TODO packed bool buffers
					global_tables_stateNow_i64[tabent.table][tabent.entry] = 1;
				}
			}
			
			// all done with message
		};
		// TODO min_delay option when no gap junctions exist
		// Also wait for recvs to finish
		// Spin it all, to probe for multimple incoming messages
		
		bool all_received = true;
		do{
			all_received = true; // at least for the empty set examined before the loop
			// TODO also try parallelizing this, perhaps?
			for( size_t idx = 0; idx < recv_off_to_node.size(); idx++ ){
				auto other_rank = recv_off_to_node.at(idx);
				const auto &recvlist_impl = engine_config.recvlist_impls.at(other_rank);
				
				if( received_sends[idx] ) continue;
				
				// otherwise it's pending
				all_received = false;
				
				
				auto &buf = recv_bufs[idx];
				auto &req = recv_requests[idx];
				
				if( received_probes[idx] ){
					// check if recv is done
					int flag = 0;
					MPI_Status status;
					MPI_Test( &req, &flag, &status);
					if( flag ){
						// received, yay !
						// Say("Recv %d.%zd", other_rank,  recvlist_impl.value_mirror_size);
						if( config.debug_netcode ){
							Say("Recv %d : %s", other_rank, NetMessage_ToString( recvlist_impl.value_mirror_size, buf).c_str());
						}
						ReceiveList( recvlist_impl, buf );
						received_sends[idx] = true;
					}
				}
				else{
					// check if probe is ready
					int flag = 0;
					MPI_Status status;
					MPI_Iprobe( other_rank, MYMPI_TAG_BUF_SEND, MPI_COMM_WORLD, &flag, &status);
					if( flag ){
						int buf_size;
						MPI_Get_count( &status, MPI_FLOAT, &buf_size );
						buf.resize( buf_size );
						PostRecv( other_rank, buf, req );
						received_probes[idx] = true;
					}
				}
				
			}
			// MPI_Finalize();
			// 	exit(1);
		} while( !all_received );
		// and clear the progress flags
		received_probes.assign( received_probes.size(), false );
		received_sends .assign( received_sends .size(), false );
		
		#endif
		
		//prepare for parallel iteration
		const float dt = engine_config.dt;
		
		// Execute all work items
		#pragma omp parallel for schedule(runtime)
		for( long long item = 0; item < engine_config.work_items; item++ ){
			if(config.debug){
				printf("item %lld start\n", item);
				// if(my_mpi.rank != 0) continue;
				// continue;
				fflush(stdout);
			}
			tabs.callbacks[item]( time, dt,
				tabs.global_constants        .data(), tabs.global_const_f32_index      [item],
				global_tables_const_f32_sizes.data(), global_tables_const_f32_arrays.data(), tabs.global_table_const_f32_index[item],
				global_tables_const_i64_sizes.data(), global_tables_const_i64_arrays.data(), tabs.global_table_const_i64_index[item],
				global_tables_state_f32_sizes.data(), global_tables_stateNow_f32, global_tables_stateNext_f32, tabs.global_table_state_f32_index[item],
				global_tables_state_i64_sizes.data(), global_tables_stateNow_i64, global_tables_stateNext_i64, tabs.global_table_state_i64_index[item],
				global_state_now, global_state_next , tabs.global_state_f32_index      [item],
				step
			);
			if(config.debug){
				printf("item %lld end\n", item);
				fflush(stdout);
			}
		}
		
		if( !initializing ){
			// output what needs to be output
			for(size_t i = 0; i < engine_config.trajectory_loggers.size(); i++){
				#ifdef USE_MPI
				assert(my_mpi.rank == 0);
				#endif
				const auto &logger = engine_config.trajectory_loggers[i];
				FILE *& fout = trajectory_open_files[i];
				
				const ScaleEntry seconds = {"sec",  0, 1.0};
				const double time_scale_factor = Scales<Time>::native.ConvertTo(1, seconds);
				
				float time_val = time * time_scale_factor;
				column_fmt.write( time_val, tmps_column );
				fprintf(fout, "%s", tmps_column);
				
				auto GetColumnValue = [
				#ifdef USE_MPI
					&engine_config,
					&global_tables_stateNow_f32,
				#endif
					&global_state_now
				]( const EngineConfig::TrajectoryLogger::LogColumn &column ){
					
					switch( column.type ){
						case EngineConfig::TrajectoryLogger::LogColumn::Type::TOPLEVEL_STATE :{
							if(column.value_type == EngineConfig::TrajectoryLogger::LogColumn::ValueType::F32){
								#ifdef USE_MPI
								if( column.on_node >= 0 && column.on_node != my_mpi.rank ){
									size_t table = engine_config.recvlist_impls.at(column.on_node).value_mirror_buffer;
									
									// scaling is done on remote node
									return global_tables_stateNow_f32[table][column.entry];
								}
								#endif
								// otherwise a local one
								return (float) (global_state_now[column.entry] * column.scaleFactor);
							}
							else if(column.value_type == EngineConfig::TrajectoryLogger::LogColumn::ValueType::I64){
								printf("but i have no flat i64 states lol\n");
								exit(2);
							}
							else{
								printf("internal error: unknown value type\n");
								exit(2);
							}
							
							break;
						}
						case EngineConfig::TrajectoryLogger::LogColumn::Type::TABLE_STATE :
						default:
							printf("internal error: unknown log type\n");
							exit(2);
					}
				};
				for( const auto &column : logger.columns ){
					char tmps_column[ column_width + 5 ];
					
					float col_val = GetColumnValue(column);
					column_fmt.write( col_val, tmps_column );
					fprintf( fout, "\t%s", tmps_column );
					// fprintf( fout, "\t%f", col_val );
					
				}
				fprintf(fout, "\n");
			}
		}
		
		// output state dump
		if(
			config.dump_raw_state_scalar
			|| config.dump_raw_state_table
		){
			if( !initializing ){
				printf("State: t = %g %s\n", time, Scales<Time>::native.name);
			}
			else{
				printf("State: t = %g %s, initialization step %lld\n", time, Scales<Time>::native.name, step);
			}
		}
		if(config.dump_raw_state_scalar){
			// print state, separated by work item
			for( size_t i = 0, itm = 1; i < state_one.size(); i++ ){
				printf("%g \t", global_state_next[i]);
				while( itm < tabs.global_state_f32_index.size() && (i + 1) == (size_t)tabs.global_state_f32_index[itm] ){
					printf("| ");
					itm++;
				}
			}
			printf("\n");
		}
		if(config.dump_raw_state_table){
			
			auto PrintVeryRawTables = [](const auto &index, const auto &arrays, const auto &sizes){
				size_t next_tabchunk = 0;
				for(size_t i = 0; i < sizes.size() ;i++){
					
					if(next_tabchunk < index.size() && i == (size_t)index.at(next_tabchunk) ){
						printf("%zd", i);
						while( next_tabchunk < index.size() && i == (size_t)index.at(next_tabchunk) ) next_tabchunk++;
					}
					printf(" \t");
					printf(" %16p \t", arrays[i]);
					for( std::size_t j = 0; j < (size_t)sizes.at(i); j++ ){
						auto val = arrays[i][j];
						printf("%s \t", (presentable_string(val)).c_str());
					}
					printf("\n");
				}
			} ;
			
			printf("RawStateF32:\n");
			PrintVeryRawTables(tabs.global_table_state_f32_index, global_tables_stateNow_f32, global_tables_state_f32_sizes);
			printf("RawStateI64:\n");
			PrintVeryRawTables(tabs.global_table_state_i64_index, global_tables_stateNow_i64, global_tables_state_i64_sizes);
			printf("RawStateNextF32:\n");
			PrintVeryRawTables(tabs.global_table_state_f32_index, global_tables_stateNext_f32, global_tables_state_f32_sizes);
			printf("RawStateNextI64:\n");
			PrintVeryRawTables(tabs.global_table_state_i64_index, global_tables_stateNext_i64, global_tables_state_i64_sizes);
		}
		
		#ifdef USE_MPI
		// wait for sends, to finish the iteration
		MPI_Waitall( send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE );
		#endif
		
		//prepare for next parallel iteration
		if( !initializing ){
			time += engine_config.dt;
		}
		
		// or a modulo-based cyclic queue LATER? will logging be so much of an issue? if so let the logger clone the state instead of duplicating the entire buffers
		std::swap(global_state_now, global_state_next); 
		std::swap(global_tables_stateNow_f32, global_tables_stateNext_f32); 
		std::swap(global_tables_stateNow_i64, global_tables_stateNext_i64); 
	}
	// done!
	
	// close loggers
	for( auto &fout : trajectory_open_files ){
		fclose(fout);
		fout = NULL;
	}
	
	gettimeofday(&run_end, NULL);
	metadata.run_time_sec = TimevalDeltaSec(run_start, run_end);
	
	
	printf("Config: %.3lf Setup: %.3lf Run: %.3lf \n", metadata.config_time_sec, metadata.init_time_sec, metadata.run_time_sec );
	#ifdef __linux__
	//get memory usage information too
	long long memResidentPeak = metadata.peak_resident_memory_bytes = getPeakResidentSetBytes();
	long long memResidentEnd = metadata.end_resident_memory_bytes = getCurrentResidentSetBytes();
	long long memHeap = getCurrentHeapBytes();
	printf("Peak: %lld Now: %lld Heap: %lld\n", memResidentPeak, memResidentEnd, memHeap );
	#endif
	//-------------------> release sim data structures, though it's not absolutely necessary at this point
	
	
	#ifdef USE_MPI
	// this is necessary, so stdio files are actually flushed
	MPI_Finalize();
	#endif
	
	return 0;
}
