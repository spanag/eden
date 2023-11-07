#include "NeuroML.h"

#include "thirdparty/pugixml-1.9/pugixml.hpp"


#include <type_traits>
#include <map>
#include <iterator>
#include <sstream>
#include <fstream>
#include <numeric>

#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <ctype.h>

//------------------> Auxiliaries

//Linear interpolation from a -> b with respective parameter from 0 -> 1
static Real Lerp(Real a, Real b, Real ratio){
	return a + (b-a)*ratio;
}

// Append vector to vector
auto AddToVec = []( auto &add_to, const auto &add_from ){
	add_to.insert( std::end(add_to), std::begin(add_from), std::end(add_from) );	
};

// Also eliminate some "../" for a better attempt to de-duplicate them.
std::string GetRelativeFilePath(const std::string &origin_path, const std::string &relative_path ){
	
	auto PathIsRelativeToCwd = []( const std::string &path ){
		#if defined _WIN32
			if( path.size() >= 1 && path[0] == '/' ) return false; // relative to current drive, or new-style absolute path
			if( path.size() >= 1
				&& (
					(path[0] >= 'A' && path[0] <= 'Z')
					|| (path[0] >= 'a' && path[0] <= 'z')
				)
				&& path[1] == ':'
			) return false; // relative to root, or cwd, of a drive
			// Does not include CP/M device names such as CON, LPT#, AUX, &c. which are erroneously considered to be files atm; Contact the devs if you need support on this
			// More details on: https://docs.microsoft.com/en-us/dotnet/standard/io/file-path-formats#identify-the-path
			return true;
		#else
			// basically Unix paths otherwise
			if( path.size() > 0 && path[0] == '/' ) return false; // is absoute path
			return true;
		#endif
	};
	
	std::string resolved_path; // the path to the file to be used, eventually
	if(!PathIsRelativeToCwd(relative_path)){
		// the path is not relative to the file including it
		resolved_path = relative_path;
	}
	else{
		// the path is relative to the directory of the file including it
		
		// generate new path, replacing current name of the file with the new href
		// contact this code's devs for more safety (not accessing sensitive files), or just run the program in chroot or something
		auto new_path = origin_path;
		// find the last slash and replace everything after it with href
		
		auto last_slash_pos = new_path.rfind('/');
		if(last_slash_pos == std::string::npos) last_slash_pos = 0;
		else last_slash_pos++; //do include the slash
		
		#if defined _WIN32
		// consider the backslash, as well
		auto last_backslash_pos = new_path.rfind('\\');
		if(last_backslash_pos == std::string::npos) last_backslash_pos = 0;
		else last_backslash_pos++; //do include the slash
		
		if( last_backslash_pos > last_slash_pos ) last_slash_pos = last_backslash_pos;
		#endif
		
		new_path.resize(last_slash_pos);
		new_path += relative_path;
		
		resolved_path = new_path;
	}
	
	// and get rid of all ../ nonsense, including symlinks? (probably not including symlinks)
	// unfortunately, realpath doesn't exist in certain environments, TODO detect and use fallback
	// const char *sResolved = realpath( new_path, NULL );
	// if( !sResolved ){
	// 	perror( ("error resolving path "+new_path).c_str() );
	// 	goto CLEANUP;
	// }
	// std::string resolved_path = sResolved;
	// free( (void *)sResolved ); sResolved = NULL;
	auto ResolvePathTextually = [  ]( const std::string &original_path ){
			
		auto tokens = string_split(std::string(original_path), "/");
		#if defined _WIN32
		// do this for Windows backslashes as well
		{
		std::vector<std::string> new_tokens;
		
		// For each token, split by slash (may contain backslashes)
		for( const auto token : tokens ){
			// For each sub-token, split by backslash
			// TODO for windows only?
			for( const auto new_token : string_split( token,  "\\") ){
				new_tokens.push_back(new_token);
			}
		}
		
		tokens = new_tokens;
		}
		#endif
		
		std::vector<std::string> resolved_tokens;
		for( int i = 0; i < (int)tokens.size(); i++ ){
			const auto &token = tokens[i];
			
			if( token == ".." ){
				if(
					!resolved_tokens.empty()
					&& *resolved_tokens.rbegin() != ".."
				){
					// move level up from last folder
					resolved_tokens.pop_back();
				}
				else{
					// higher than cwd, probably
					resolved_tokens.push_back(token);
				}
			}
			else if( token == "." ){
				// skip
			}
			else if( token == "" ){
				if( i == 0 ){
					// path begins with a slash, keep for the sake of Unix
					resolved_tokens.push_back(token);
				}
				else{
					// superfluous slashes, skip
				}
			}
			else{
				resolved_tokens.push_back(token);
			}
		}
		
		std::string resolved_path;						
		for( int i = 0; i < (int)resolved_tokens.size(); i++ ){
			const auto &token = resolved_tokens[i];
			resolved_path += token;
			if( i + 1 < (int)resolved_tokens.size() ) resolved_path += "/";
		}
		return resolved_path;
	};
	resolved_path = ResolvePathTextually( resolved_path );

	return resolved_path;
}

//------------------> XML Utilities start

void ReportErrorInFile(FILE *error_log, const char *filename, const pugi::xml_node &node, const char *format, ...){
	va_list args;
	va_start(args, format);
	ReportErrorInFile_Base(error_log, filename, node.offset_debug(), format, args);
	va_end (args);
}

// The internal state of XML files being opened to perform an import. Is required by the import process logger.
struct NmlFileContext{
	std::string filename;
	
	pugi::xml_document doc;
	// TODO add a offset <-> source line cache, to not scan the file over and over
};
struct NmlImportContext{
	std::list<NmlFileContext> documents_opened;
	
	// keyed by internal pointer of each document's root element
	// (the api was made so elegant that this is the only way to retrieve associate xml_node with element)
	std::unordered_map<const void *, const NmlFileContext *> files_by_root_element; 
	
	//NmlImportContext(){}
};

//The object handling smart logging, with references to the specific point in the specific file being referred to in the log entry
struct ImportLogger{
	
	const NmlImportContext &import_context;
	FILE *error_log;
	
	// control this one from general logging, why not?
	bool debug_printing;
	
	// keep a line position cache for multiple warnings, instead of bailing out LATER
	
	ImportLogger(const NmlImportContext &_import, FILE *_error_log = stderr, bool _debug = false)
	: import_context(_import), error_log(_error_log), debug_printing(_debug) {  }
	
	// NOTE since only the DOM element pointer can determine the opened file, the DOM trees should be kept in ptr-preserving containers (like std::list)
	
	// get filename from XML element, NULL if missing/unknown
	const char *GetFilenameFromElement( const pugi::xml_node &node ) const {
		const char *filename = NULL;
		const void *context_key = node.root().internal_object();
		//printf("Retrieve key %p\n",context_key);
		if(import_context.files_by_root_element.count(context_key)){
			filename = import_context.files_by_root_element.at(context_key)->filename.c_str();
		}
		return filename;		
	}
	
	// LATER maybe add RAII context stack (e.g. in file, in segment, ...)
	void _doit(const char *filename, const ptrdiff_t file_byte_offset, const char *format, va_list args) const {
		ReportErrorInFile_Base(error_log, filename, file_byte_offset, format, args);
	}
	void _doit(const pugi::xml_node &node, const char *format, va_list args) const {
		// ReportErrorInFile_Base(error_log, GetFilenameFromElement(node), node.offset_debug(), format, args);
		_doit(GetFilenameFromElement(node), node.offset_debug(), format, args);
	}
	
	void error(const char *filename, const ptrdiff_t file_byte_offset, const char *format, ...) const {
		va_list args;
		va_start(args, format);
		_doit(filename, file_byte_offset, format, args);
		va_end (args);
	}
	void error(const pugi::xml_node &node, const char *format, ...) const {
		va_list args;
		va_start(args, format);
		
		_doit(node, format, args);
		
		va_end (args);
	}
	void error(va_list args, const char *filename, const ptrdiff_t file_byte_offset, const char *format) const {
		_doit(filename, file_byte_offset, format, args);
	}
	//same as error, really
	void warning(const char *filename, const ptrdiff_t file_byte_offset, const char *format, ...) const {
		va_list args;
		va_start(args, format);
		_doit(filename, file_byte_offset, format, args);
		va_end (args);
	}
	void warning(const pugi::xml_node &node, const char *format, ...) const {
		va_list args;
		va_start(args, format);
		
		_doit(node, format, args);
		
		va_end (args);
	}
	void warning(va_list args, const char *filename, const ptrdiff_t file_byte_offset, const char *format) const {
		_doit(filename, file_byte_offset, format, args);
	}
	void warn_experimental( const pugi::xml_node &node, const char *featname ) const {
		warning(node, "%s is an experimental extension; its behaviour may change in the future.", featname); // TODO add hashset of such things to whine once for each type.
	}
};

struct LogWithElement : public ILogProxy{
	const ImportLogger &log;
	const pugi::xml_node &eLement;
	LogWithElement(const ImportLogger &_l, const pugi::xml_node &_e): log(_l), eLement(_e) {}
	void error(const char *format, ...) const {
		va_list args;
		va_start(args, format);
		log._doit(eLement, format, args);
		va_end (args);
	}
	void warning(const char *format, ...) const {
		va_list args;
		va_start(args, format);
		log._doit(eLement, format, args);
		va_end (args);
	}
};

//------------------> XML Utilities end

//------------------> Dimensions

const char * Dimensionless::NAME = "Dimensionless";
const char * Time::NAME = "Time";
const char * Frequency::NAME = "Frequency";
const char * Length::NAME = "Length";
const char * Voltage::NAME = "Voltage";
const char * Current::NAME = "Current";
const char * Resistance::NAME = "Resistance";
const char * Conductance::NAME = "Conductance";
const char * Resistivity::NAME = "Resistivity";
const char * Conductivity::NAME = "Conductivity";
const char * Capacitance::NAME = "Capacitance";
const char * SpecificCapacitance::NAME = "Specific Capacitance";
const char * Concentration::NAME = "Concentration";
const char * Permeability::NAME = "Permeability";
const char * RhoFactor::NAME = "Rho Factor";
const char * Temperature::NAME = "Temperature";

//------------------> Scales of Dimensions

template<> const ScaleEntry Scales<Dimensionless>::native = {"units", 0, 1.0};

template<> const ScaleEntry Scales<Time>::native = {"us", -6, 1.0};
template<> const ScaleList Scales<Time>::scales = {
	{"hour", 3600, 1.0}, // some funny people really use this
	{"s"  ,  0, 1.0},
	{"ms" , -3, 1.0},
	{"msec", -3, 1.0}, // a cute alias
	{"us" , -6, 1.0},
	{"usec", -6, 1.0}, // a cute alias
	{"sec",0, 1.0}, // why was this not used
	{"second",0, 1.0},
	{"seconds",0, 1.0}
};
// TODO validate unit scale conversions in LEMS !!!
template<> const ScaleEntry Scales<Frequency>::native = {"MHz", +6, 1.0};
template<> const ScaleList Scales<Frequency>::scales = {
	{"per_s"  ,  0, 1.0},
	{"Hz"     ,  0, 1.0},
	{"per_ms" , +3, 1.0},
	{"kHz"    , +3, 1.0},
	{"per_us" , +6, 1.0},
	{"MHz"    , +6, 1.0},
};

template<> const ScaleEntry Scales<Length>::native = {"um", -6, 1.0};
template<> const ScaleList Scales<Length>::scales = {
	{"m" ,  0, 1.0},
	{"cm", -2, 1.0}, // physiological
	{"mm", -3, 1.0},
	{"um", -6, 1.0},
};

template<> const ScaleEntry Scales<Voltage>::native = {"mV", -3, 1.0};
template<> const ScaleList Scales<Voltage>::scales = {
	{"V" ,  0, 1.0},
	{"mV", -3, 1.0},
};

template<> const ScaleEntry Scales<Current>::native = {"nA", -9, 1.0};
template<> const ScaleList Scales<Current>::scales = {
	{"A" ,   0, 1.0},
	{"mA", - 3, 1.0},
	{"uA", - 6, 1.0},
	{"nA", - 9, 1.0},
	{"pA", -12, 1.0},
};

template<> const ScaleEntry Scales<Resistance>::native = {"Mohm", 6, 1.0};
template<> const ScaleList Scales<Resistance>::scales = {
	{"ohm" , 0, 1.0},
	{"kohm", 3, 1.0}, // physioogical
	{"Mohm", 6, 1.0},
};

template<> const ScaleEntry Scales<Conductance>::native = {"uS", -6, 1.0};
template<> const ScaleList Scales<Conductance>::scales = {
	{"S" ,   0, 1.0},
	{"mS", - 3, 1.0},
	{"uS", - 6, 1.0},
	{"nS", - 9, 1.0},
	{"pS", -12, 1.0},
};

template<> const ScaleEntry Scales<Resistivity>::native = {"kohm_cm", 1, 1.0};
template<> const ScaleList Scales<Resistivity>::scales = {
	{"kohm_cm" ,   1, 1.0},
	{"ohm_m"   ,   0, 1.0},
	{"ohm_cm"  , - 2, 1.0},
};

template<> const ScaleEntry Scales<Conductivity>::native = {"mS_per_cm2", 1, 1.0};
template<> const ScaleList Scales<Conductivity>::scales = {
	{"S_per_m2"  , 0, 1.0},
	{"mS_per_cm2", 1, 1.0},
	{"S_per_cm2" , 4, 1.0},
	{"S_per_mm2" , 6, 1.0},
};

template<> const ScaleEntry Scales<Capacitance>::native = {"pF",-12, 1.0};
template<> const ScaleList Scales<Capacitance>::scales = {
	{"F" ,  0, 1.0},
	{"mF",- 3, 1.0},
	{"uF",- 6, 1.0},
	{"nF",- 9, 1.0},
	{"pF",-12, 1.0},
};

template<> const ScaleEntry Scales<SpecificCapacitance>::native = {"uF_per_cm2", -2, 1.0};
template<> const ScaleList Scales<SpecificCapacitance>::scales = {
	{"uF_per_cm2" ,-2, 1.0},
	{"F_per_m2"   , 0, 1.0},
	{"uF_per_mm2" , 0, 1.0},
};

template<> const ScaleEntry Scales<Concentration>::native = {"uM" , -3, 1.0};
template<> const ScaleList Scales<Concentration>::scales = {
	{"mol_per_m3" ,  0, 1.0},
	{"mol_per_cm3",  6, 1.0},
	{"M"          ,  3, 1.0},
	{"mM"         ,  0, 1.0},
	{"uM"         , -3, 1.0},
	{"nM"         , -6, 1.0},
};

template<> const ScaleEntry Scales<Permeability>::native = {"m_per_s", 0, 1.0};
template<> const ScaleList Scales<Permeability>::scales = {
	{"m_per_s"  , 0, 1.0},
	{"cm_per_s" ,-2, 1.0},
	{"um_per_ms",-3, 1.0},
	{"cm_per_ms", 1, 1.0},
};

template<> const ScaleEntry Scales<RhoFactor>::native = {"mol_per_m_per_A_per_s" , 0, 1.0};
template<> const ScaleList Scales<RhoFactor>::scales = {
	{"mol_per_m_per_A_per_s"   ,  0, 1.0},
	{"mol_per_cm_per_uA_per_ms", 11, 1.0},
};

template<> const ScaleEntry Scales<Temperature>::native = {"K" , 0, 1.0};
template<> const ScaleList Scales<Temperature>::scales = {
	{"K"   ,  0, 1.0},
	{"degK",  0, 1.0},
	{"degC",  0, 1.0, 273.15},
	{"degF",  0, 5.0/9.0,  ( -32 * 5.0/9.0 ) + 273.15},
};

//                                          M   L   T   I   K   N   J
const Dimension LEMS_Pure               = { 0,  0,  0,  0,  0,  0,  0};
const Dimension LEMS_Time               = { 0,  0, +1,  0,  0,  0,  0};
const Dimension LEMS_Frequency          = { 0,  0, -1,  0,  0,  0,  0};
const Dimension LEMS_Voltage            = { 1,  2, -3, -1,  0,  0,  0};
const Dimension LEMS_Area               = { 0,  2,  0,  0,  0,  0,  0};
const Dimension LEMS_Concentration      = { 0, -3,  0,  0,  0,  1,  0};
const Dimension LEMS_Current            = { 0,  0,  0,  1,  0,  0,  0};
const Dimension LEMS_CurrentDensity     = { 0, -2,  0,  1,  0,  0,  0};
const Dimension LEMS_Temperature        = { 0,  0,  0,  0,  1,  0,  0};
const Dimension LEMS_Conductance        = {-1, -2,  3,  2,  0,  0,  0};
const Dimension LEMS_Conductivity       = {-1, -4,  3,  2,  0,  0,  0};
const Dimension LEMS_Capacitance        = {-1, -2,  4,  2,  0,  0,  0};
/* To obtain a consistent units system, solve the following linear equation between exponents:

	  M  L  T  I  K  N  J
	/ 0  0  1  0  0  0  0 \   / native Mass      (-18) \   / preferred Time          (usec      ) \
	| 1  2 -3 -1  0  0  0 |   | native Length    (- 6) |   | preferred Voltage       (millivolts) |
	| 0  0  0  1  0  0  0 |   | native Time      (- 6) |   | preferred Current       (nanoamps  ) |
	| 0  1  0  0  0  0  0 | * | native Current   (- 9) | = | preferred Length        (microns   ) |
	| 0  0  0  0  1  0  0 |   | native Kelvins   (  0) |   | preferred Kelvins       (kelvins   ) |
	| 0 -3  0  0  0  1  0 |   | native Substance (- 3) |   | preferred Concentration (micromolar) |
	\ 0  0  0  0  0  0  1 /   \ native J         (  0) /   \ preferred J             (who knows ) /	
*/
void DimensionSet::AddDefaults(){
	//             M  L  T  I  K  N  J
	Add(Dimension( 0, 0, 0, 0, 0, 0, 0), {"unitless"               , { "units",   0}, {{ "units",   0}} });
	Add(Dimension( 0, 0, 0, 0, 0, 0, 0), {"none"                   , });
	Add(Dimension( 0, 0, 1, 0, 0, 0, 0), DimensionInfo::FromOld<Time               >("time"                   ));
	Add(Dimension( 0, 0,-1, 0, 0, 0, 0), DimensionInfo::FromOld<Frequency          >("per_time"               ));
	Add(Dimension( 1, 2,-3,-1, 0, 0, 0), DimensionInfo::FromOld<Voltage            >("voltage"                ));
	Add(Dimension(-1,-2, 3, 2, 0, 0, 0), DimensionInfo::FromOld<Conductance        >("conductance"            ));
	Add(Dimension(-1,-4, 3, 2, 0, 0, 0), DimensionInfo::FromOld<Conductivity       >("conductanceDensity"     ));
	Add(Dimension(-1,-2, 4, 2, 0, 0, 0), DimensionInfo::FromOld<Capacitance        >("capacitance"            ));
	Add(Dimension(-1,-4, 4, 2, 0, 0, 0), DimensionInfo::FromOld<SpecificCapacitance>("specificCapacitance"    ));
	Add(Dimension( 1, 2,-3,-2, 0, 0, 0), DimensionInfo::FromOld<Resistance         >("resistance"             ));
	Add(Dimension( 1, 3,-3,-2, 0, 0, 0), DimensionInfo::FromOld<Resistivity        >("resistivity"            ));
	Add(Dimension( 0, 0, 0, 1, 0, 0, 0), DimensionInfo::FromOld<Current            >("current"                ));
	Add(Dimension( 0, 1, 0, 0, 0, 0, 0), DimensionInfo::FromOld<Length             >("length"                 ));
	Add(Dimension( 0,-3, 0, 0, 0, 1, 0), DimensionInfo::FromOld<Concentration      >("concentration"          ));
	Add(Dimension( 0, 1,-1, 0, 0, 0, 0), DimensionInfo::FromOld<Permeability       >("permeability"           ));
	Add(Dimension( 0, 0, 0, 0, 1, 0, 0), DimensionInfo::FromOld<Temperature        >("temperature"            ));
	Add(Dimension( 0,-1,-1,-1, 0, 1, 0), DimensionInfo::FromOld<RhoFactor          >("rho_factor"             ));
	
	Add(Dimension(-1,-3, 3, 1, 0, 0, 0), {"per_voltage"            , { "per_mV"         ,   3 }, { { "per_V"          , 0 }, { "per_mV"        ,   3 } } });
	Add(Dimension( 0, 0, 1, 1, 0, 0, 0), {"charge"                 , { "pC"             , -12 }, { { "C"              , 0 }, { "pC"            , -12 } } });
	Add(Dimension( 0, 0, 1, 1, 0,-1, 0), {"charge_per_mole"        , { "C_per_mol"      ,   0 }, { { "C_per_mol"      , 0 }, { "nA_ms_per_amol",   6 } } });
	Add(Dimension( 0,-2, 0, 1, 0, 0, 0), {"currentDensity"         , { "uA_per_cm2"     , - 2 }, { { "A_per_m2"       , 0 }, { "uA_per_cm2"    , - 2 }, { "mA_per_cm2",   1 } } });
	Add(Dimension( 0, 2, 0, 0, 0, 0, 0), {"area"                   , { "um2"            , -12 }, { { "m2"             , 0 }, { "cm2"           , - 4 }, { "um2"       , -12 } } });
	Add(Dimension( 0, 3, 0, 0, 0, 0, 0), {"volume"                 , { "um3"            , -18 }, { { "m3"             , 0 }, { "cm3"           , - 6 }, { "litre"     , - 3 }, { "um3", -18 } } });
	Add(Dimension( 0, 0, 0, 0, 0, 1, 0), {"substance"              , { "amol"           , -18 }, { { "mol"            , 0 }, { "amol"          , -18 } } });
	Add(Dimension( 1, 2,-2, 0,-1,-1, 0), {"idealGasConstantDims"   , { "J_per_K_per_mol",   0 }, { { "J_per_K_per_mol", 0 }  } });
	Add(Dimension(-2,-4, 6, 3, 0, 0, 0), {"conductance_per_voltage", { "nS_per_mV"      , - 6 }, { { "S_per_V"        , 0 }, { "nS_per_mV"     , - 6 } } });
	
	Add(Dimension( 1, 2,-4,-1, 0, 0, 0), {"voltage_per_time"       , { "mV_per_usec"      ,   3 }, { { "V_per_sec"      , 0 }, { "mV_per_msec"   ,   0 }, { "mV_per_usec"   ,   3 } } });	
}

//------------------> Misc utilities

// map core components to LEMS core implementations
const char *LEMS_CoreComponents_filename = "LEMS_CoreComponents.inc.xml";


//------------------> Model object code

bool ComponentType::GetExposureAndDimension( Int exp_seq, Dimension &dim_out ) const {
	if( exp_seq >= 0 ){
		dim_out = getExposureDimension( exp_seq );
		return true;
	}
	else return false;
}
bool ComponentType::GetCurrentOutputAndDimension( Dimension &dim_out ) const {
	return GetExposureAndDimension( common_exposures.current, dim_out );
}
bool ComponentType::GetVoltageRequirementAndDimension( Dimension &dim_out ) const {
	return GetRequirementAndDimension( common_requirements.membrane_voltage, dim_out );
}

bool InputSource::HasSpikeOut( const CollectionWithNames<ComponentType> &component_types ) const {
	if( type == COMPONENT ){
		if( component_types.get( component.id_seq ).common_event_outputs.spike_out >= 0 ) return true;
		else return false;
	}
	else{
		if(
			type == SPIKE_LIST ||  type == SPIKE_PERIODIC || type == SPIKE_RANDOM || type == SPIKE_POISSON || type == SPIKE_POISSON_REF || type == PYNN_SPIKE_POISSON
			|| type == TIMED_SYNAPTIC || type == POISSON_SYNAPSE || type == POISSON_SYNAPSE_TRANSIENT
		) return true;
		// TODO check composite ones
		else return false;
	}
}
bool InputSource::GetCurrentOutputAndDimension( const CollectionWithNames<ComponentType> &component_types, const CollectionWithNames<SynapticComponent> &synapses, Dimension &dim_out ) const {
	if( type == COMPONENT ){
		return component_types.get( component.id_seq ).GetCurrentOutputAndDimension( dim_out );
	}
	else{
		if(
			type == PULSE || type == SINE ||  type == RAMP || type == VOLTAGE_CLAMP || type == VOLTAGE_CLAMP_TRIPLE
			|| type == TIMED_SYNAPTIC || type == POISSON_SYNAPSE || type == POISSON_SYNAPSE_TRANSIENT || type == COMPOSITE
		){
			dim_out = LEMS_Current;
			return true;
		}
		else if(
			type == PULSE_DL || type == SINE_DL ||  type == RAMP_DL || type == COMPOSITE_DL
		){
			dim_out = Dimension::Unity();
			return true;
		}
		else return false;
	}
}
bool InputSource::GetVoltageInputAndDimension( const CollectionWithNames<ComponentType> &component_types, const CollectionWithNames<SynapticComponent> &synapses, Dimension &dim_out ) const {
	if( type == COMPONENT ){
		return component_types.get( component.id_seq ).GetVoltageRequirementAndDimension( dim_out );
	}
	else{
		if( type == VOLTAGE_CLAMP ||  type == VOLTAGE_CLAMP_TRIPLE ){
			dim_out = LEMS_Voltage;
			return true;
		}
		else if( type == TIMED_SYNAPTIC || type == POISSON_SYNAPSE || type == POISSON_SYNAPSE_TRANSIENT ){
			return synapses.get(synapse).GetVoltageInputAndDimension( component_types, dim_out );
		}
		else return false;
		// TODO composite inputs
	}
}

bool SynapticComponent::HasSpikeIn( const CollectionWithNames<ComponentType> &component_types ) const {
	if( type == COMPONENT ){
		if( component_types.get( component.id_seq ).common_event_inputs.spike_in >= 0 ) return true;
		else return false;
	}
	else{
		if(
			type == ALPHA || type == EXP || type == EXPTWO || type == EXPTHREE || type == ALPHA_CURRENT || type == BLOCKING_PLASTIC
			|| type == PYNN_EXP_COND || type == PYNN_ALPHA_COND || type == PYNN_EXP_CURR || type == PYNN_ALPHA_CURR
		) return true;
		else return false;
		// TODO composite synapses
	}
}
bool SynapticComponent::GetCurrentOutputAndDimension( const CollectionWithNames<ComponentType> &component_types, Dimension &dim_out ) const {
	if( type == COMPONENT ){
		return component_types.get( component.id_seq ).GetCurrentOutputAndDimension( dim_out );
	}
	else{
		if(
			type == SILENT
			|| type == ALPHA ||  type == EXP || type == EXPTWO || type == EXPTHREE || type == ALPHA_CURRENT || type == BLOCKING_PLASTIC
			|| type == GAP ||  type == GAPHALF || type == GRADED_PRINZ_ET_AL
			|| type == PYNN_EXP_COND || type == PYNN_ALPHA_COND || type == PYNN_EXP_CURR || type == PYNN_ALPHA_CURR
		){
			dim_out = LEMS_Current;
			return true;
		}
		else return false;
	}
	// TODO composite synapses
}
bool SynapticComponent::GetVoltageInputAndDimension( const CollectionWithNames<ComponentType> &component_types, Dimension &dim_out ) const {
	if( type == COMPONENT ){
		return component_types.get( component.id_seq ).GetVoltageRequirementAndDimension( dim_out );
	}
	else{
		if(
			type == SILENT
			|| type == ALPHA ||  type == EXP || type == EXPTWO || type == EXPTHREE || type == ALPHA_CURRENT || type == BLOCKING_PLASTIC
			|| type == GAP ||  type == GAPHALF || type == GRADED_PRINZ_ET_AL
			|| type == PYNN_EXP_COND || type == PYNN_ALPHA_COND || type == PYNN_EXP_CURR || type == PYNN_ALPHA_CURR
		){
			dim_out = LEMS_Voltage;
			return true;
		}
		else return false;
	}
}
// should be removed because we always care about voltage dimensionality? no, not necessarily
bool SynapticComponent::HasVpeer( const CollectionWithNames<ComponentType> &component_types ) const {
	Dimension dummy;
	return GetVpeerInputAndDimension(component_types, dummy);
}
bool SynapticComponent::GetVpeerInputAndDimension( const CollectionWithNames<ComponentType> &component_types, Dimension &dim_out ) const {
	if( type == COMPONENT ){
		const auto &comptype = component_types.get( component.id_seq );
		return comptype.GetRequirementAndDimension( comptype.common_requirements.peer_voltage, dim_out );
	}
	else{
		if(
			type == SILENT
			|| type == GAP ||  type == GAPHALF || type == GRADED_PRINZ_ET_AL
		){
			dim_out = LEMS_Voltage;
			return true;
		}
		else return false;
	}
}

bool ArtificialCell::GetCurrentInputAndDimension( const CollectionWithNames<ComponentType> &component_types, Dimension &dim_out ) const {
	if( type == COMPONENT ){
		const auto &comptype = component_types.get( component.id_seq );
		return comptype.GetRequirementAndDimension( comptype.common_requirements.external_current, dim_out );
	}
	else if( type == SPIKE_SOURCE ){
		return false; // perhaps something else LATER
	}
	else{
		if(
			type == IAF || type == IAF_REF || type == IZH_2007 || type == ADEX
			|| type == PYNN_IF_CURR_ALPHA || type == PYNN_IF_CURR_EXP
			|| type == PYNN_IF_COND_ALPHA || type == PYNN_IF_COND_EXP
			|| type == PYNN_EIF_COND_EXP_ISFA_ISTA || type == PYNN_EIF_COND_ALPHA_ISFA_ISTA
			|| type == PYNN_HH_COND_EXP
		){
			dim_out = LEMS_Current;
			return true;
		}
		else if(
			type == IZH
		){
			dim_out = Dimension::Unity();
			return true;
		}
		else if(
			type == IAF_TAU || type == IAF_TAU_REF || type == FN || type == FN_1969 || type == PINSKY_RINZEL_CA3
		){
			return false;
		}
		else return false;
	}
}
bool ArtificialCell::GetVoltageExposureAndDimension( const CollectionWithNames<ComponentType> &component_types, Dimension &dim_out ) const {
	if( type == COMPONENT ){
		const auto &comptype = component_types.get( component.id_seq );
		return comptype.GetExposureAndDimension( comptype.common_exposures.membrane_voltage, dim_out );
	}
	else if( type == SPIKE_SOURCE ){
		return false; // perhaps something else LATER
	}
	else{
		if(
			type == IAF || type == IAF_REF || type == IZH_2007 || type == ADEX || type == IZH
			|| type == IAF_TAU || type == IAF_TAU_REF || type == PINSKY_RINZEL_CA3
			|| type == PYNN_IF_CURR_ALPHA || type == PYNN_IF_CURR_EXP
			|| type == PYNN_IF_COND_ALPHA || type == PYNN_IF_COND_EXP
			|| type == PYNN_EIF_COND_EXP_ISFA_ISTA || type == PYNN_EIF_COND_ALPHA_ISFA_ISTA
			|| type == PYNN_HH_COND_EXP
			
		){
			dim_out = LEMS_Voltage;
			return true;
		}
		else if(
			type == FN || type == FN_1969
		){
			dim_out = Dimension::Unity();
			return true;
		}
		else return false; // TODO DL core components
	}
}

bool ArtificialCell::HasSpikeIn( const CollectionWithNames<ComponentType> &component_types ) const {
	if( type == COMPONENT ){
		if( component_types.get( component.id_seq ).common_event_inputs.spike_in >= 0 ) return true;
		else return false;
	}
	else if( type == SPIKE_SOURCE ){
		return false; // perhaps something else LATER
	}
	else{
		if(
			type == IAF || type == IAF_REF || type == IAF_TAU || type == IAF_TAU_REF
			|| type == IZH || type == IZH_2007 || type == ADEX
			|| type == PYNN_IF_CURR_ALPHA || type == PYNN_IF_CURR_EXP
			|| type == PYNN_IF_COND_ALPHA || type == PYNN_IF_COND_EXP
			|| type == PYNN_EIF_COND_EXP_ISFA_ISTA || type == PYNN_EIF_COND_ALPHA_ISFA_ISTA
			|| type == FN || type == FN_1969 || type == PINSKY_RINZEL_CA3 || type == PYNN_HH_COND_EXP
		) return false; // no core type
		else return false;
	}
}
bool ArtificialCell::HasSpikeOut( const CollectionWithNames<ComponentType> &component_types, const CollectionWithNames<InputSource> &input_sources ) const {
	if( type == COMPONENT ){
		if( component_types.get( component.id_seq ).common_event_outputs.spike_out >= 0 ) return true;
		else return false;
	}
	else if( type == SPIKE_SOURCE ){
		return input_sources.get(spike_source_seq).HasSpikeOut(component_types);
	}
	else{
		if(
			type == IAF || type == IAF_REF || type == IAF_TAU || type == IAF_TAU_REF
			|| type == IZH || type == IZH_2007 || type == ADEX
			|| type == PYNN_IF_CURR_ALPHA || type == PYNN_IF_CURR_EXP
			|| type == PYNN_IF_COND_ALPHA || type == PYNN_IF_COND_EXP
			|| type == PYNN_EIF_COND_EXP_ISFA_ISTA || type == PYNN_EIF_COND_ALPHA_ISFA_ISTA
		) return true;
		else if(
			type == FN || type == FN_1969 || type == PINSKY_RINZEL_CA3 || type == PYNN_HH_COND_EXP
		) return false;
		else return false;
	}
}

bool CellType::GetCurrentInputAndDimension( const CollectionWithNames<ComponentType> &component_types, Dimension &dim_out ) const {
	if( type == ARTIFICIAL ){
		return artificial.GetCurrentInputAndDimension( component_types, dim_out );
	}
	else{
		// it's a physical cell
		dim_out = LEMS_Current;
		return true;
	}
}
bool CellType::GetVoltageExposureAndDimension( const CollectionWithNames<ComponentType> &component_types, Dimension &dim_out ) const {
	if( type == ARTIFICIAL ){
		return artificial.GetVoltageExposureAndDimension( component_types, dim_out );
	}
	else{
		// it's a physical cell
		dim_out = LEMS_Voltage;
		return true;
	}
}
bool CellType::HasSpikeIn( const CollectionWithNames<ComponentType> &component_types ) const {
	if( type == ARTIFICIAL ){
		return artificial.HasSpikeIn( component_types );
	}
	else{
		// physical cells do not receive spikes
		return false;
	}
}
bool CellType::HasSpikeOut( const CollectionWithNames<ComponentType> &component_types, const CollectionWithNames<InputSource> &input_sources ) const {
	if( type == ARTIFICIAL ){
		return artificial.HasSpikeOut( component_types, input_sources );
	}
	else{
		// physical cells emit spikes, in principle; check if particulat compartment has defined Vt further on
		return true;
	}
}

const NameMap< Int ComponentType::CommonRequirements::* > ComponentType::CommonRequirements::names = {
	{"time"				, &ComponentType::CommonRequirements::time }, 
	{"temperature"		, &ComponentType::CommonRequirements::temperature }, 
	{"v"				, &ComponentType::CommonRequirements::membrane_voltage }, // typically physical
	{"V"				, &ComponentType::CommonRequirements::membrane_voltage }, // typically dimensionless
	{"surfaceArea"		, &ComponentType::CommonRequirements::membrane_surface_area }, 
	{"iSyn"		        , &ComponentType::CommonRequirements::external_current }, // typically physical
	{"ISyn"		        , &ComponentType::CommonRequirements::external_current }, // typically dimensionless
	{"iCa"				, &ComponentType::CommonRequirements::ion_current },
	{"iCa2"				, &ComponentType::CommonRequirements::ion_current }, // XXX handle it needing both at the same time
	{"initialConcentration", &ComponentType::CommonRequirements::initial_concentration_intra }, 
	{"initialExtConcentration", &ComponentType::CommonRequirements::initial_concentration_extra }, 
	{"caConc"			, &ComponentType::CommonRequirements::calcium_concentration_intra }, 
	{"alpha"			, &ComponentType::CommonRequirements::gate_rate_alpha }, 
	{"beta"				, &ComponentType::CommonRequirements::gate_rate_beta  }, 
	{"rateScale"		, &ComponentType::CommonRequirements::gate_rate_scale }, 
	{"vpeer"			, &ComponentType::CommonRequirements::peer_voltage },
	{"blockFactor"	, &ComponentType::CommonRequirements::block_factor },
	{"plasticityFactor"	, &ComponentType::CommonRequirements::plasticity_factor },
};
const NameMap< Int ComponentType::CommonExposures::* > ComponentType::CommonExposures::names = {
	{"v"				, &ComponentType::CommonExposures::membrane_voltage }, // typically physical
	{"V"				, &ComponentType::CommonExposures::membrane_voltage }, // typically dimensionless
	{"i"				, &ComponentType::CommonExposures::current }, // typically physical
	{"I"				, &ComponentType::CommonExposures::current }, // typically dimensionless
	{"g"				, &ComponentType::CommonExposures::conductance }, // typically dimensionless
	{"r"				, &ComponentType::CommonExposures::rate }, 
	{"t"				, &ComponentType::CommonExposures::time }, 
	{"x"				, &ComponentType::CommonExposures::variable }, 
	{"q"				, &ComponentType::CommonExposures::gate_variable }, 
	{"fcond"			, &ComponentType::CommonExposures::conductivity_fraction }, // when it's an ion channel gate
	{"fopen"			, &ComponentType::CommonExposures::conductivity_fraction },	// when it's the whole ion channel
	{"concentration"    , &ComponentType::CommonExposures::concentration_intra },
	{"extConcentration"	, &ComponentType::CommonExposures::concentration_extra },
	{"factor"			, &ComponentType::CommonExposures::scaling_factor },
	{"blockFactor"		, &ComponentType::CommonExposures::block_factor },
	{"plasticityFactor"	, &ComponentType::CommonExposures::plasticity_factor },
};
const NameMap< Int ComponentType::CommonEventInputs::* > ComponentType::CommonEventInputs::names = {
	{"in"    , &ComponentType::CommonEventInputs::spike_in },
};
const NameMap< Int ComponentType::CommonEventOutputs::* > ComponentType::CommonEventOutputs::names = {
	{"spike" , &ComponentType::CommonEventOutputs::spike_out },
};

// due to silly reasons an implementation of even constexpr static member variables must be declared in a cpp file to link, fixed in c++17
constexpr int Dimension::* Dimension::members[];
constexpr Int ComponentType::CommonRequirements::* ComponentType::CommonRequirements::members[];
constexpr Int ComponentType::CommonExposures   ::* ComponentType::CommonExposures   ::members[];
constexpr Int ComponentType::CommonEventInputs ::* ComponentType::CommonEventInputs ::members[];
constexpr Int ComponentType::CommonEventOutputs::* ComponentType::CommonEventOutputs::members[];

bool ParseSynapseLocation(const ILogProxy &log, Network::Projection::Connection::Type type, const char *sLocation,  Simulation::LemsQuantityPath::SynapsePath::Location &loca){
	typedef Simulation::LemsQuantityPath::SynapsePath::Location SynLoc;
	typedef Network::Projection::Connection::Type ConnType;
	if(     strcmp(sLocation, "post") == 0) loca = SynLoc::POST;
	else if(strcmp(sLocation, "pre" ) == 0) loca = SynLoc::PRE;
	else{
		log.error("unknown synapse location %s", sLocation);
		return false;
	}
	// check if location is appropriate
	if( (loca == SynLoc::PRE) && !(type == ConnType::ELECTRICAL || type == ConnType::CONTINUOUS) ){
		log.error("'pre' specifier can only be used for graded synapses", sLocation);
		return false;
	}
	if( (loca == SynLoc::POST) && !(type == ConnType::ELECTRICAL || type == ConnType::CONTINUOUS || type == ConnType::SPIKING) ){
		log.error("'post' specifier can only be used for event-based or graded synapses", sLocation);
		return false;
	}
	return true;
}
void Model::LemsSegmentLocatorToString( const Network &net, const Simulation::LemsSegmentLocator &path, std::string &ret) const {
	const Network::Population &population = net.populations.get(path.population);
	const CellType &cell = cell_types.get(population.component_cell);
	ret += net.populations.getName(path.population) + ("["+accurate_string(path.cell_instance)+"]");
	// if it's a physical cell, add an explicit segment locator, just to be nice
	if(cell.type == CellType::PHYSICAL){
		// add segment locator
		ret += "/";
		ret += accurate_string(path.segment_seq);
		if(path.fractionAlong != 0.5) ret += accurate_string(path.fractionAlong);
		// don't add a slash in the end, refer to the very point on cell
	}
	return;
}

// helpers for LEMS id's
bool Model::GetLemsQuantityPathType_FromLems(const Simulation::LemsInstanceQuantityPath &path, const ComponentInstance &compinst, ComponentType::NamespaceThing::Type &type, Dimension &dimension) const {
	if(!compinst.ok()) return false;
	Int comp_type_seq = compinst.id_seq;
	
	const ComponentType &comp_type = component_types.get(comp_type_seq);
	const auto &namespace_thing_seq = path.namespace_thing_seq;
	type = comp_type.name_space.get(namespace_thing_seq).type;
	dimension = comp_type.getNamespaceEntryDimension(namespace_thing_seq);
	return true;
}
bool Model::GetLemsQuantityPathType_SynapticComponent(const Simulation::SynapticComponentQuantityPath &path, const SynapticComponent &syn, ComponentType::NamespaceThing::Type &type, Dimension &dimension) const {
	typedef SynapticComponent Syn;
	typedef Simulation::SynapticComponentQuantityPath SynPath;
	
	if(path.type == SynPath::LEMS){
		return GetLemsQuantityPathType_FromLems(path.lems_quantity_path, syn.component, type, dimension);
	}
	else if(path.type == SynPath::BLOCK){
		if(syn.blopla.block_mechanism.type == SynapticComponent::BlockingPlasticSynapse::BlockMechanism::NONE){ assert(false); return false; }
		// NB they are all lemsified for now
		return GetLemsQuantityPathType_FromLems(path.block.lems_quantity_path, syn.blopla.block_mechanism.component, type, dimension);
	}
	else if(path.type == SynPath::PLASTICITY){
		if(syn.blopla.plasticity_mechanism.type == SynapticComponent::BlockingPlasticSynapse::PlasticityMechanism::NONE){ assert(false); return false; }
		// NB they are all lemsified for now
		return GetLemsQuantityPathType_FromLems(path.plasticity.lems_quantity_path, syn.blopla.plasticity_mechanism.component, type, dimension);
	}
	else if(path.type == SynPath::NATIVE){
		// single exp is not lemsified
		if( syn.type == Syn::EXP || syn.type == Syn::GAP ){
			if( path.native_entry == SynPath::GBASE ){
				type = ComponentType::NamespaceThing::PROPERTY;
				dimension = LEMS_Conductance;
				return true;
			}
			else if( path.native_entry == SynPath::G ){
				type = ComponentType::NamespaceThing::STATE;
				dimension = LEMS_Conductance;
				return true;
			}
			else if( path.native_entry == SynPath::EREV ){
				type = ComponentType::NamespaceThing::PROPERTY;
				dimension = LEMS_Voltage;
				return true;
			}
			else if( path.native_entry == SynPath::TAU ){
				type = ComponentType::NamespaceThing::PROPERTY;
				dimension = LEMS_Time;
				return true;
			}
			else{ assert(false); return false; }
		}
		// LATER double synapses and more!
		else{ assert(false); return false; }
	}
	else{ assert(false); return false; }
}
bool Model::GetLemsQuantityPathType_InputInstance(const Simulation::InputInstanceQuantityPath &path, const InputSource &input, ComponentType::NamespaceThing::Type &type, Dimension &dimension) const {
	typedef Simulation::InputInstanceQuantityPath InpPath;
	// TODO work with child synapses in input component!
	if(path.type == InpPath::LEMS){
		return GetLemsQuantityPathType_FromLems(path.lems_quantity_path, input.component, type, dimension);
	}
	else if(path.type == InpPath::SYNAPSE){
		return GetLemsQuantityPathType_SynapticComponent(path.synapse, synaptic_components.get(input.synapse), type,dimension);
	}
	else if(path.type == InpPath::NATIVE){
		// dc pulse is not lemsified
		if(input.type == InputSource::Type::PULSE){
			if( path.native_entry == InpPath::AMPLITUDE ){
				type = ComponentType::NamespaceThing::PROPERTY;
				dimension = LEMS_Current;
				return true;
			}
			else if( path.native_entry == InpPath::DURATION
				|| path.native_entry == InpPath::DELAY ){
				type = ComponentType::NamespaceThing::PROPERTY;
				dimension = LEMS_Time;
				return true;
			}
			else{ assert(false); return false; }
		}
		// FIXME inputs with synapse component!
		// FIXME spike sources!
		else{ assert(false); return false; }
	}
	else{ assert(false); return false; }
}
bool Model::GetLemsQuantityPathType(const Network &net, const Simulation::LemsQuantityPath &path, ComponentType::NamespaceThing::Type &type, Dimension &dimension) const {
	typedef ComponentType::NamespaceThing NamespaceThing;
	
	if(path.type == Simulation::LemsQuantityPath::CELL){
		const Network::Population &population = net.populations.get(path.population);
		const CellType &_cell_type = cell_types.get(population.component_cell);
		
		if(_cell_type.type != CellType::ARTIFICIAL){assert(false); return false;} // otherwise it would be a segment path
		const auto &cell = _cell_type.artificial;
		
		const auto &path_cell = path.cell;
		
		if(path_cell.type == Simulation::LemsQuantityPath::CellPath::INPUT){
			const auto &input = input_sources.get( cell.spike_source_seq );
			return GetLemsQuantityPathType_InputInstance(path_cell.input, input, type, dimension);
		}
		else{
			// NB all artificial cells should be lemsified atm
			return GetLemsQuantityPathType_FromLems(path_cell.lems_quantity_path, cell.component, type, dimension);
		}
	}
	else if(path.type == Simulation::LemsQuantityPath::SEGMENT){
		const Network::Population &population = net.populations.get(path.population);
		const CellType &_cell_type = cell_types.get(population.component_cell);
		
		assert(_cell_type.type == CellType::PHYSICAL);
		const auto &cell = _cell_type.physical; (void) cell; // TODO
		
		typedef Simulation::LemsQuantityPath::SegmentPath SegmentPath;
		switch(path.segment.type){
			case SegmentPath::VOLTAGE       : type = NamespaceThing::STATE; dimension = LEMS_Voltage      ; return true;
			case SegmentPath::CALCIUM_INTRA : type = NamespaceThing::STATE; dimension = LEMS_Concentration; return true; // TODO check again it it's correct
			case SegmentPath::CALCIUM2_INTRA: type = NamespaceThing::STATE; dimension = LEMS_Concentration; return true;
			case SegmentPath::CALCIUM_EXTRA : type = NamespaceThing::STATE; dimension = LEMS_Concentration; return true; // FIXME
			case SegmentPath::CALCIUM2_EXTRA: type = NamespaceThing::STATE; dimension = LEMS_Concentration; return true;
			default: assert(false); return false;
		}
	}
	else if(path.type == Simulation::LemsQuantityPath::CHANNEL){
		const Network::Population &population = net.populations.get(path.population);
		const CellType &_cell_type = cell_types.get(population.component_cell);
		assert(_cell_type.type == CellType::PHYSICAL);
		const PhysicalCell &cell = _cell_type.physical;
		
		// const Morphology &morph = morphologies.get(cell_type.physical.morphology);
		const BiophysicalProperties &bioph = biophysics.at(cell.biophysicalProperties);
		
		const auto &dist = bioph.membraneProperties.channel_specs.at(path.channel.distribution_seq);
		const auto &channel = ion_channels.get(dist.ion_channel);
		// lots TODO here
		
		typedef Simulation::LemsQuantityPath::ChannelPath ChannelPath;
		switch(path.channel.type){
			// these are, affected by gbar iirc (TODO check if a lems channel can give pure current, probably)
			case ChannelPath::I        : type = NamespaceThing::DERIVED; dimension = LEMS_Current     ; return true;
			case ChannelPath::G        : type = NamespaceThing::DERIVED; dimension = LEMS_Conductance; return true; // TODO check again if it's correct
			case ChannelPath::I_DENSITY: type = NamespaceThing::DERIVED; dimension = LEMS_CurrentDensity; return true;
			case ChannelPath::G_DENSITY: type = NamespaceThing::DERIVED; dimension = LEMS_Conductivity; return true;
			// these are about the channel mechanism
			case ChannelPath::Q        : type = NamespaceThing::STATE  ; dimension = LEMS_Pure; return true; // FIXME that's not always the case!
			case ChannelPath::SUBQ     : type = NamespaceThing::STATE  ; dimension = LEMS_Pure; return true; // FIXME that's not always the case!
			case ChannelPath::LEMS     : return GetLemsQuantityPathType_FromLems( path.channel.lems_quantity_path, channel.component, type, dimension);
			default: assert(false); return false;
		}
	}/*
	else if(path.type == Simulation::LemsQuantityPath::ION_POOL){
		
	}*/
	else if(path.type == Simulation::LemsQuantityPath::SYNAPSE){
		const auto &path_synapse = path.synapse;
		if(path_synapse.proj_seq < 0) return false;
		// NB connections in a projection are not guaranteed to use the same component;
		// LATER constrain them to guarantee this.
		if(path_synapse.conn_seq < 0) return false;
		const Network::Projection &proj = net.projections.get(path_synapse.proj_seq);
		const Network::Projection::Connection &conn = proj.connections.atSeq(path_synapse.conn_seq);
		Int syn_seq = path_synapse.GetSynSeq(conn);
		const SynapticComponent &syn = synaptic_components.get(syn_seq);
		
		return GetLemsQuantityPathType_SynapticComponent(path_synapse, syn, type, dimension);
	}
	else if(path.type == Simulation::LemsQuantityPath::INPUT){
		const Network::InputList &list = net.input_lists.get(path.input.list_seq);
		const InputSource &input = input_sources.get(list.component);
		
		return GetLemsQuantityPathType_InputInstance(path.input, input, type, dimension);
	}
	// experimental extensions
	else if(path.type == Simulation::LemsQuantityPath::DATAREADER){
		const Network::TimeSeriesReader &reader = net.data_readers.get(path.reader.read_seq);
		const Network::TimeSeriesReader::InputColumn &column = reader.columns.get(path.reader.colu_seq);
		type = NamespaceThing::STATE; dimension = column.dimension;
		return true;
	}
	else if(path.type == Simulation::LemsQuantityPath::NETWORK){
		// not supported yet!
		return false;
	}
	else{
		printf("path type %d\n", (int) path.type);
		assert(false);
		return false;
	}
}
bool Model::GetLemsEventPathType(const Simulation::LemsInstanceEventPath &path, const ComponentInstance &compinst, Simulation::LemsInstanceEventPath::Type &type) const {
	type = path.type;
	return true;
}
bool Model::GetLemsEventPathType(const Simulation::InputInstanceEventPath &path, const InputSource &input, Simulation::LemsInstanceEventPath::Type &type) const {
	typedef Simulation::InputInstanceEventPath InpPath;
	// TODO work with child synapses in input component!
	if(path.type == InpPath::LEMS){
		return GetLemsEventPathType(path.lems_event_path, input.component, type);
	}
	else if(path.type == InpPath::NATIVE){
		// spikeArray is not lemsified
		if( path.native_entry == InpPath::SPIKE ){
			type = Simulation::LemsInstanceEventPath::OUT;
			return true;
		}
		else{ assert(false); return false; }
	}
	else{ assert(false); return false; }
}
bool Model::GetLemsEventPathType(const Network &net, const Simulation::LemsEventPath &path, Simulation::LemsInstanceEventPath::Type &type) const {
	typedef Simulation::LemsInstanceEventPath::Type Type;
	
	if(path.type == Simulation::LemsEventPath::CELL){
		const Network::Population &population = net.populations.get(path.population);
		const CellType &_cell_type = cell_types.get(population.component_cell);
		
		if(_cell_type.type != CellType::ARTIFICIAL){assert(false); return false;} // otherwise it would be a segment path
		const auto &cell = _cell_type.artificial;
		
		const auto &path_cell = path.cell;
		
		if(path_cell.type == Simulation::LemsEventPath::Cell::INPUT){
			const auto &input = input_sources.get( cell.spike_source_seq );
			return GetLemsEventPathType(path_cell.input, input, type);
		}
		else{
			// NB all artificial cells should be lemsified atm
			return GetLemsEventPathType(path_cell.lems_event_path, cell.component, type);
		}
	}
	else if(path.type == Simulation::LemsEventPath::SEGMENT){
		const Network::Population &population = net.populations.get(path.population);
		const CellType &_cell_type = cell_types.get(population.component_cell);
		
		assert(_cell_type.type == CellType::PHYSICAL);
		const auto &cell = _cell_type.physical; (void) cell; // TODO
		
		typedef Simulation::LemsEventPath::Segment SegmentPath;
		switch(path.segment.type){
			case SegmentPath::SPIKE       : type = Type::OUT; return true;
			default: assert(false); return false;
		}
	}
	// experimental extensions
	else if(path.type == Simulation::LemsEventPath::EVENTREADER){
		type = Type::OUT;
		return true;
	}
	else{
		printf("event path type %d\n", (int) path.type);
		assert(false);
		return false;
	}
}
bool Model::LemsQuantityPathToString(const ComponentInstance &inst, const Simulation::LemsInstanceQuantityPath &path, std::string &ret) const {
	assert(inst.ok());
	Int comptype_seq = inst.id_seq;
	assert(component_types.has(comptype_seq));
	const auto &comptype = component_types.get(comptype_seq);
	
	const char *sName = comptype.name_space.getName(path.namespace_thing_seq);
	if(!sName) return false;
	ret += sName;
	// all namespace things are leaves, until LATER
	return true;
}
bool Model::LemsQuantityPathToString(const SynapticComponent &syn, const Simulation::SynapticComponentQuantityPath &path, std::string &ret) const  {
	typedef Simulation::SynapticComponentQuantityPath Path;
	if( path.type == Path::LEMS ){
		return LemsQuantityPathToString(syn.component, path.lems_quantity_path, ret);
	}
	else if( path.type == Path::NATIVE ){
		typedef Path::NativeEntry Nat;
		switch( path.native_entry ){
			case Nat::GBASE : ret += "gbase"    ; return true;
			case Nat::EREV  : ret += "erev"     ; return true;
			case Nat::TAU   : ret += "tauDecay" ; return true;
			case Nat::G     : ret += "g"        ; return true;
			default: return false;
		}
	}
	else if( path.type == Path::BLOCK ){
		ret += "blockMechanisms[0]/";
		if(path.block.type == Path::Block::LEMS){
			return LemsQuantityPathToString(syn.blopla.block_mechanism.component, path.block.lems_quantity_path, ret);
		}
		else return false;
	}
	else if( path.type == Path::PLASTICITY ){
		ret += "plasticityMechanisms[0]/";
		if(path.plasticity.type == Path::Plasticity::LEMS){
			return LemsQuantityPathToString(syn.blopla.plasticity_mechanism.component, path.plasticity.lems_quantity_path, ret);
		}
		else return false;
	}
	else return false;
}
bool Model::LemsQuantityPathToString(const InputSource &input, const Simulation::InputInstanceQuantityPath &path, std::string &ret) const {
	typedef Simulation::InputInstanceQuantityPath Path;
	if( path.type == Path::LEMS ){
		return LemsQuantityPathToString(input.component, path.lems_quantity_path, ret);
	}
	else if( path.type == Path::NATIVE ){
		switch( path.native_entry ){
			case Path::AMPLITUDE : ret += "amplitude"; return true;
			case Path::DURATION  : ret += "duration" ; return true;
			case Path::DELAY     : ret += "delay"    ; return true;
			default: return false;
		}
	}
	else if( path.type == Path::SYNAPSE ){
		ret += "synapse/";
		return LemsQuantityPathToString(synaptic_components.get(input.synapse), path.synapse, ret);
	}
	else return false;
}
bool Model::LemsQuantityPathToString(const ArtificialCell &cell, const Simulation::LemsQuantityPath::CellPath &path, std::string &ret) const {
	typedef Simulation::LemsQuantityPath::CellPath Path;
	if( path.type == Path::LEMS ){
		return LemsQuantityPathToString(cell.component, path.lems_quantity_path, ret);
	}
	else if( path.type == Path::INPUT ){
		assert( cell.type == ArtificialCell::SPIKE_SOURCE );
		assert( cell.spike_source_seq >= 0 );
		const auto &input = input_sources.get( cell.spike_source_seq );
		return LemsQuantityPathToString(input, path.input, ret);
	}
	else return false;
}
bool Model::LemsQuantityPathToString(const Network &net, const Simulation::LemsQuantityPath &path, std::string &ret) const {
	const auto &type = path.type;
	typedef Simulation::LemsQuantityPath Path;
	
	ret.clear();
	if(type == Path::CELL
	|| type == Path::SEGMENT
	|| type == Path::CHANNEL
	|| type == Path::ION_POOL
	){
		const Network::Population &population = net.populations.get(path.population);
		const CellType &cell = cell_types.get(population.component_cell);
		LemsSegmentLocatorToString(net, path, ret); ret += "/";
		
		if(type == Path::CELL){
			assert(cell.type == CellType::ARTIFICIAL);
			return LemsQuantityPathToString(cell.artificial, path.cell, ret);
		}
		else if(type == Path::SEGMENT){
			const auto type = path.segment.type;
			typedef Path::SegmentPath Path;
			if     ( type == Path::VOLTAGE ) ret += "v";
			else if( type == Path::CALCIUM_INTRA ) ret += "caConc";
			else if( type == Path::CALCIUM2_INTRA ) ret += "caConc2";
			else return false;
			return true;
		}
		else if(type == Path::CHANNEL){
			ret += "biophysicalProperties/membraneProperties/";
			const auto &bioph = biophysics.at(cell.physical.biophysicalProperties);
			const auto &chandist = bioph.membraneProperties.channel_specs[path.channel.distribution_seq];
			const auto &channel = ion_channels.get(chandist.ion_channel);
			
			ret += chandist.name; ret += "/";
			typedef Path::ChannelPath Path; 
			switch(path.channel.type){
				case Path::I : ret += "i"; return true;
				case Path::G : ret += "g"; return true;
				case Path::I_DENSITY : ret += "iDensity"; return true;
				case Path::G_DENSITY : ret += "gDensity"; return true;
				case Path::Q :{
					ret += ion_channels.getName(chandist.ion_channel); ret += "/";
					ret += channel.gates.getName(path.channel.gate_seq); ret += "/q";
					return true;
				}
				// TODO SUBQ, LEMS
				default: break;
			}
			return false; // none of supported path types, somehow
		}
		// TODO ionpool!
		else return false;
	}
	else if(type == Path::SYNAPSE){
		ret += net.projections.getName(path.synapse.proj_seq) + ("["+accurate_string(path.synapse.conn_seq)+"]/");
		const Network::Projection &proj = net.projections.get(path.synapse.proj_seq);
		const Network::Projection::Connection &conn = proj.connections.contents[path.synapse.conn_seq];
		
		if(path.synapse.location == Path::SynapsePath::Location::POST) ret += "post/";
		else if(path.synapse.location == Path::SynapsePath::Location::PRE) ret += "pre/";
		else return false;
		
		Int syn_seq = path.synapse.GetSynSeq(conn);
		const SynapticComponent &syn = synaptic_components.get(syn_seq);
		// LATER access weight, delay as well?
			
		return LemsQuantityPathToString(syn, path.synapse, ret);
	}
	else if(type == Path::INPUT){
		ret += net.input_lists.getName(path.input.list_seq) + ("["+accurate_string(path.input.inst_seq)+"]/");
		const Network::InputList &list = net.input_lists.get(path.input.list_seq);
		const InputSource &input = input_sources.get(list.component);
		return LemsQuantityPathToString(input, path.input, ret);
	}
	// experimental extensions
	else if(type == Path::DATAREADER){
		ret += net.data_readers.getName(path.reader.read_seq) + ("["+accurate_string(path.reader.inst_seq)+"]/") + net.data_readers.get(path.reader.read_seq).columns.getName(path.reader.colu_seq);
		return true;
	}
	else{
		printf("path to string: type %d not supported yet\n", (int)type);
		return false;
	}
}
bool Model::LemsEventPathToString(const ComponentInstance &inst, const Simulation::LemsInstanceEventPath &path, std::string &ret) const {
	typedef Simulation::LemsInstanceEventPath Path;
	assert(inst.ok());
	Int comptype_seq = inst.id_seq;
	assert(component_types.has(comptype_seq));
	const auto &comptype = component_types.get(comptype_seq);
	
	if(path.type == Path::Type::IN ) ret += comptype.event_inputs .getName(path.event_port_seq);
	else if(path.type == Path::Type::OUT) ret += comptype.event_outputs.getName(path.event_port_seq);
	else{ assert(false); return false;}
	// all namespace things are leaves, until LATER
	return true;
}
bool Model::LemsEventPathToString(const InputSource &input, const Simulation::InputInstanceEventPath &path, std::string &ret) const {
	typedef Simulation::InputInstanceEventPath Path;
	if( path.type == Path::LEMS ){
		return LemsEventPathToString(input.component, path.lems_event_path, ret);
	}
	else if( path.type == Path::NATIVE ){
		switch( path.native_entry ){
			case Path::SPIKE : ret += "spike"; return true;
			default: return false;
		}
	}
	else return false;
}
bool Model::LemsEventPathToString(const ArtificialCell &cell, const Simulation::LemsEventPath::Cell &path, std::string &ret) const {
	typedef Simulation::LemsEventPath::Cell Path;
	if( path.type == Path::LEMS ){
		return LemsEventPathToString(cell.component, path.lems_event_path, ret);
	}
	else if( path.type == Path::INPUT ){
		assert( cell.type == ArtificialCell::SPIKE_SOURCE );
		assert( cell.spike_source_seq >= 0 );
		const auto &input = input_sources.get( cell.spike_source_seq );
		return LemsEventPathToString(input, path.input, ret);
	}
	else return false;
}
bool Model::LemsEventPathToString(const Network &net, const Simulation::LemsEventPath &path, std::string &ret) const {
	const auto &type = path.type;
	typedef Simulation::LemsEventPath Path;
	
	ret.clear();
	if(type == Path::CELL
	|| type == Path::SEGMENT
	){
		const Network::Population &population = net.populations.get(path.population);
		const CellType &cell = cell_types.get(population.component_cell);
		LemsSegmentLocatorToString(net, path, ret); ret += "/";
		
		if(type == Path::CELL){
			assert(cell.type == CellType::ARTIFICIAL);
			return LemsEventPathToString(cell.artificial, path.cell, ret);
		}
		else if(type == Path::SEGMENT){
			const auto type = path.segment.type;
			typedef Path::Segment Path;
			if     ( type == Path::SPIKE ) ret += "spike";
			else return false;
			return true;
		}
		else return false;
	}
	// experimental extensions
	else if(type == Path::EVENTREADER){
		ret += net.event_readers.getName(path.reader.read_seq) + ("["+accurate_string(path.reader.inst_seq)+"]/") + net.event_readers.get(path.reader.read_seq).ports.getName(path.reader.port_seq);
		return true;
	}
	else{
		printf("path to string: type %d not supported yet\n", (int)type);
		return false;
	}
}
// parse according to the definitive specification: the source code
// https://github.com/NeuroML/jNeuroML/issues/50
// https://github.com/NeuroML/org.neuroml.export/blob/master/src/main/java/org/neuroml/export/utils/LEMSQuantityPath.java#L136
// https://github.com/NeuroML/org.neuroml.export/blob/master/src/main/java/org/neuroml/export/neuron/LEMSQuantityPathNeuron.java#L118
// https://github.com/NeuroML/org.neuroml.export/blob/master/src/main/java/org/neuroml/export/neuron/NeuronWriter.java#L1110
// https://github.com/NeuroML/org.neuroml.export/blob/development/src/main/java/org/neuroml/export/neuron/NeuronWriter.java#L654
// note that LEMS does fully modular resolution (but does not support multiple compartments) :
// https://github.com/LEMS/jLEMS/blob/development/src/main/java/org/lemsml/jlems/core/sim/RunnableAccessor.java#L29

// <population name>[instance]/[<redundant population>/]v for point neurons
// <population name>/<instance>/<cell type, redundant>/<segment>/v
// ... <segment>/<biophysical properties name>/membraneProperties/<ion channel distribution>/[i, g, iDensity, gDensity]
// ... <segment>/<biophysical properties name>/membraneProperties/<ion channel distribution>/<channel type>/<gate>/q
// <postsynaptic population name>/<instance>/<cell type, redundant>/<segment>/synapses:<synapse type>:<incrementing instance of same type on same segment of same cell>/g 
// 	see https://github.com/NeuroML/org.neuroml.export/blob/master/src/main/java/org/neuroml/export/neuron/NeuronWriter.java#L581
// template<auto whatever> typename IdToSeq;
template<typename V, Int (V::*prop) >
Int IdToSeq(const CollectionWithNames<V> &container, Int id){
	return container.*prop.getSequential(id);
};
// template<auto whatever> struct SeqToSeq;
template<typename V>
struct SeqToSeq{
	using Prop = Int V::*;
	const Prop prop;
	SeqToSeq(const Prop _p):prop(_p){}
	Int operator()(const V &ob, Int id) const { 
		return ( 0 <= id && id < ob.*prop ) ? id : -1 ;
	}
};
// template<typename V, typename S, S (V::*prop) >
// template<typename V, typename S>
// template<S (V::*prop)>
// Int SeqToSeq(const V &collection, Int id){
// 	return ( 0 <= id && id < (Int)(collection.*prop) ) ? id : -1 ;
// };
template<typename GroupCollection, typename Get_Seq>
bool ParseLemsGroupLocator(const ILogProxy &log, const std::vector<std::string> tokens, const char *group_type_name, const GroupCollection &group_collection, const Get_Seq &GetSeq, Int &group_seq, Int &instance_seq, Int &tokens_consumed) {
	// note: tokens are spearted by '/' in the following
	Int pop_token_id = tokens_consumed;
	if(tokens_consumed >= (int) tokens.size()){
		return false;
	}
	const std::string &pop_token = tokens[pop_token_id];
	std::string instance_string;
	
	size_t bracketIndex =  pop_token.find_first_of("[");
	std::string pop_name = pop_token.substr(0,bracketIndex);
	
	group_seq = group_collection.get_id(pop_name.c_str());
	if(group_seq < 0){
		log.error("target %s %s not found", group_type_name, pop_name.c_str());
		return false;
	}
	const auto &group = group_collection.get(group_seq);

	if(pop_name.length() == pop_token.length()){
		// no bracket, use pop/id/... instead
		if( pop_token_id + 1 >= (int)tokens.size() ){
			log.error("not enough factors for instance ID");
			return false;
		}
		instance_string = tokens[pop_token_id + 1];
		tokens_consumed = pop_token_id + 2;
	}
	else{
		//with bracket, pop[id]/...
		instance_string = pop_token.substr(bracketIndex + 1, pop_token.find_first_of("]") - (bracketIndex + 1));
		tokens_consumed = pop_token_id + 1;
	}
	
	Int instance_id;
	if( !StrToL(instance_string.c_str(), instance_id) ){
		log.error("target instance \"%s\" not an integer number", instance_string.c_str());
		return false;
	}
	instance_seq = GetSeq(group, instance_id);
	if(instance_seq < 0){
		log.error("target instance %s not found in %s", instance_string.c_str(), pop_name.c_str());
		return false;
	}
	return true;
}
bool ParseLemsCellLocator(const ILogProxy &log, const std::vector<std::string> tokens, const Network &net, Simulation::LemsSegmentLocator &path, Int &tokens_consumed) {
	return ParseLemsGroupLocator(log, tokens, "population", net.populations, [](const Network::Population &group, Int id){return group.instances.getSequential(id);}, path.population, path.cell_instance, tokens_consumed);
}
bool Model::ParseLemsSegmentLocator(const ILogProxy &log, const std::vector<std::string> tokens, const Network &net, Simulation::LemsSegmentLocator &path, Int &tokens_consumed) const {
	
	if( tokens.size() < 1 ){
		log.error("target path must have at least 1 slash-delimited factor");
		return false;
	}
	if(!( ParseLemsCellLocator(log, tokens, net, path, tokens_consumed) )) return false;
	
	// TODO pop might not be set when selecting populations or projections !
	const Network::Population &population = net.populations.get(path.population);
	
	int segment_token_id = tokens_consumed;
	if(segment_token_id < (int)tokens.size()){
		
		const std::string &segment_token = tokens[segment_token_id];
		// check if the segment token is further below because this is the redundant pop name
		if( segment_token == cell_types.getName(population.component_cell) ){
			//printf("skip redundant %s\n", segment_token.c_str());
			segment_token_id++;
		}
		
	}
	int segprop_token_id = -1;
	auto MatchesFractionalSegmentId = [](const char *sId){
		int dots_counted = 0;
		for(const char *s = sId; *s; s++){
			char c = *s;
			if(isdigit(c)) ; // continue
			else if(c == '.') dots_counted++;
			else return false;
		}
		return dots_counted <= 1;
	};
	auto GetFractionalSegmentId = [&MatchesFractionalSegmentId](const char *sId, Int &segment_id, Real &fractionAlong){
		if(!MatchesFractionalSegmentId(sId)) return false;
		auto dot = strchr(sId,'.');
		if(dot){
			if(!StrToL(sId, segment_id, false)) return false; // allow stopping at dot
			if(!( StrToF(dot, fractionAlong) && 0 <= fractionAlong && fractionAlong <= 1 )) return false;
		}
		else{
			if(!StrToL(sId, segment_id)) return false;
			fractionAlong = 0.5; // when not mentioned
		}
		return true;
	};
	Int segment_id = -1;
	// if(!( segment_token_id < (int)tokens.size() && StrToL(tokens[segment_token_id].c_str(), segment_id) )){
	if(!( segment_token_id < (int)tokens.size() && GetFractionalSegmentId(tokens[segment_token_id].c_str(), segment_id, path.fractionAlong) )){
		// could be a single compartment cell, or an artificial cell
		segment_id = -1; //set it to missing
		// this token is a property of the single compartment cell, ignore
		segprop_token_id = segment_token_id;
	}
	else{
		// accept token
		segprop_token_id = segment_token_id + 1;
	}
	
	// check if segment ID exists, or can be defaulted if unset
	const CellType &cell_type = cell_types.get(population.component_cell);
	if( cell_type.type == CellType::ARTIFICIAL ){
		if( segment_id < 0 ){
			// set it to default
			segment_id = 0;
		}
		if( segment_id != 0 ){
			log.error("artificial cell only has segment 0");
			return false;
		}
		
		path.segment_seq = 0;
		path.fractionAlong = 0.5; // for good measure
	}
	else if(cell_type.type == CellType::PHYSICAL ){
		
		const Morphology &morph = morphologies.get(cell_type.physical.morphology);
		if( segment_id < 0 ){
			if( morph.segments.size() == 1 ){
				// set it to default
				segment_id = 0;
				path.fractionAlong = 0.5;
			}
			else{
				std::string complaint = "target path needs segment ID, because cell has multiple segments. Setting to implicit default: segment ID = 0";
				// add diagnostic, may be useful? 
				if( 0 <= segprop_token_id && segprop_token_id < (int)tokens.size() ) complaint += "(path following cell id was \""+tokens[segment_token_id]+"\")";
				log.warning("%s", complaint.c_str());
				// return false;
				// TODO perhaps stop complaining after warning too many times
				
				segment_id = 0;
				path.fractionAlong = 0.5;
			}
		}
		
		path.segment_seq = morph.segments.getSequential(segment_id);
		if(path.segment_seq < 0){
			log.error("target segment %s not found in cell", tokens[segment_token_id].c_str());
			return false;
		}
	}
	else{
		log.error("internal error: LEMS segment locator: cell type type %d", cell_type.type);
		return false;
	}
	
	tokens_consumed = segprop_token_id;
	return true;
}
bool Model::ParseLemsQuantityPathInComponent(const ILogProxy &log, const ComponentInstance &instance, const std::vector<std::string> &tokens, Simulation::LemsInstanceQuantityPath &lems_instance_qty_path, Int &tokens_consumed ) const {
	
	if( (int)tokens.size() <= tokens_consumed ){
		log.error("path needs to specify a property of LEMS component %s", ( tokens.empty() ? "" : tokens[tokens.size() - 1].c_str() ) );
		return false;
	}
	if( tokens_consumed + 1 < (int)tokens.size() ){
		log.error("LEMS child component quantities not yet supported" );
		return false;
	}
	const char *propname = tokens[tokens_consumed].c_str();
	tokens_consumed++;
	
	// get comp.type
	Int comptype_seq = instance.id_seq;
	if( !component_types.has(comptype_seq) ){
		log.error("internal error: LEMS quantity path missing component type %d", (int) instance.id_seq);
		return false;
	}
	
	const auto &comptype = component_types.get(comptype_seq);
	
	Int &namespace_id = lems_instance_qty_path.namespace_thing_seq;
	namespace_id = comptype.name_space.get_id(propname);
	if( namespace_id < 0 ){
		log.error("%s is not a define d quantity in component type %s", propname, component_types.getName(comptype_seq));
		return false;
	}
	
	return true;
}
bool Model::ParseLemsQuantityPath_SynapticComponent(const ILogProxy &log, const SynapticComponent &syn, const std::vector<std::string> &tokens, Simulation::SynapticComponentQuantityPath &path, Int &tokens_consumed ) const {
	typedef SynapticComponent Syn;
	typedef Simulation::SynapticComponentQuantityPath SynPath;
	int mainprop_token_id = tokens_consumed;
	
	// now branch according to property type
	if(mainprop_token_id >= (Int)tokens.size()){
		log.error("not enough factors for synaptic component property");
		return false;
	}
	
	const std::string &mainprop = tokens[mainprop_token_id];
	
	// TODO weight?
	// work with LEMSified stuff as a priority
	if( syn.type != Syn::BLOCKING_PLASTIC && syn.component.ok() ){
		path.type = SynPath::LEMS;
		if( !ParseLemsQuantityPathInComponent(log, syn.component, tokens, path.lems_quantity_path, tokens_consumed ) ) return false;
		return true;
	}
	else if( syn.type == Syn::BLOCKING_PLASTIC ){
		// LATER children will be supported...
		// TODO can they be a multitude?
		if(mainprop == "blockMechanisms[0]"){
			typedef Syn::BlockingPlasticSynapse::BlockMechanism Blo;
			const Blo &blo = syn.blopla.block_mechanism;
			if(blo.type == Blo::Type::NONE){
				log.error("%s is not present on this component", mainprop.c_str());
				return false;
			}
			tokens_consumed++;
			path.type = SynPath::BLOCK;
			
			if(!(tokens_consumed+1 <= (Int)tokens.size())){
				log.error("incomplete path for %s", mainprop.c_str());
				return false;
			}
			
			assert(blo.component.ok());
			path.block.type = SynPath::Block::LEMS;
			return ParseLemsQuantityPathInComponent(log, blo.component, tokens, path.block.lems_quantity_path, tokens_consumed );
		}
		else if(mainprop == "plasticityMechanisms[0]"){
			typedef Syn::BlockingPlasticSynapse::PlasticityMechanism Pla;
			const Pla &pla = syn.blopla.plasticity_mechanism;
			if(pla.type == Pla::Type::NONE){
				log.error("%s is not present on this component", mainprop.c_str());
				return false;
			}
			tokens_consumed++;
			path.type = SynPath::PLASTICITY;
			
			if(!(tokens_consumed+1 <= (Int)tokens.size())){
				log.error("incomplete path for %s", mainprop.c_str());
				return false;
			}
			
			assert(pla.component.ok());
			path.plasticity.type = SynPath::Plasticity::LEMS;
			return ParseLemsQuantityPathInComponent(log, pla.component, tokens, path.plasticity.lems_quantity_path, tokens_consumed );
		}
		else{
			// check for the main, synaptic component which happens to be lemsified
			return ParseLemsQuantityPathInComponent(log, syn.component, tokens, path.lems_quantity_path, tokens_consumed );
		} 
	}
	
	// otherwise it's a native, or complicated, type.
	path.type = SynPath::NATIVE;
	
	tokens_consumed++; // note: this token will be used in the following, unshift or refactor if necessary later
	auto &native_entry = path.native_entry;
	typedef SynPath::NativeEntry Nat;
	
	// is it a terminal property? (to avoid repeated checks in the following)
	if(tokens_consumed == (Int)tokens.size()){
		if(     mainprop == "gbase"      ){native_entry = Nat::GBASE; return true; }
		else if(mainprop == "conductance"){native_entry = Nat::GBASE; return true; }
		else if(mainprop == "erev"       ){native_entry = Nat::EREV ; return true; }
		else if(mainprop == "tauDecay"   ){native_entry = Nat::TAU  ; return true; }
		else if(mainprop == "g"          ){native_entry = Nat::G    ; return true; }
		else{
			log.error("unknown leaf property %s", mainprop.c_str()); // LATER log in token possibly in attribute
			return false;
		}
	}
	else{
		log.error("unknown non-leaf property %s", mainprop.c_str());
		return false;
	}
	
	if( syn.type == Syn::EXP || syn.type == Syn::GAP ){
		if(( syn.type == Syn::EXP && !( native_entry == Nat::GBASE || native_entry == Nat::EREV || native_entry == Nat::TAU || native_entry == Nat::G ) )
		|| ( syn.type == Syn::GAP && !( native_entry == Nat::GBASE ) )
		){
			log.error("property %s does not match synaptic component type", mainprop.c_str()); return false;
		}
		return true; // all checks passed
	}
	// LATER composite?
	else{
		log.error("synaptic component type not supported yet");
		return false;
	}
}
bool Model::ParseLemsQuantityPath_InputInstance(const ILogProxy &log, const InputSource &input, const std::vector<std::string> &tokens, Simulation::InputInstanceQuantityPath &path, Int &tokens_consumed ) const {
	typedef InputSource::Type InpType;
	typedef Simulation::InputInstanceQuantityPath InpPath;
	// TODO weight?
	int mainprop_token_id = tokens_consumed;
	
	// now branch according to property type
	if(mainprop_token_id >= (Int)tokens.size()){
		log.error("not enough factors for input property");
		return false;
	}
	
	const std::string &mainprop = tokens[mainprop_token_id];
	
	// const bool is_native = ( input.type == InpType::PULSE );
	const bool has_synapse = (
		input.type == InpType::TIMED_SYNAPTIC
		|| input.type == InpType::POISSON_SYNAPSE
		|| input.type == InpType::POISSON_SYNAPSE_TRANSIENT
	);
	
	// synapse containing stuff are more than just lemsified
	if( has_synapse && mainprop == "synapse" ){
		tokens_consumed++;
		path.type = InpPath::SYNAPSE;
		return ParseLemsQuantityPath_SynapticComponent(log, synaptic_components.get(input.synapse), tokens, path.synapse, tokens_consumed);
	}
	
	// work with simple LEMSified stuff as a priority
	if( input.component.ok() ){
		path.type = InpPath::LEMS;
		if( !ParseLemsQuantityPathInComponent(log, input.component, tokens, path.lems_quantity_path, tokens_consumed ) ) return false;
		return true;
	}
	
	// otherwise it's a native, or complicated, type.
	path.type = InpPath::NATIVE;
	tokens_consumed++; // note: this token will be used in the following, unshift or refactor if necessary later
	
	if(input.type == InpType::PULSE){
		// is it a terminal property? (to avoid repeated checks in the following)
		if(tokens_consumed == (Int)tokens.size()){
			if(     mainprop == "amplitude"){ path.native_entry = InpPath::AMPLITUDE; return true; }
			else if(mainprop == "duration" ){ path.native_entry = InpPath::DURATION ; return true; }
			else if(mainprop == "delay"    ){ path.native_entry = InpPath::DELAY    ; return true; }
			else{
				log.error("unknown leaf property %s", mainprop.c_str()); // LATER log in token possibly in attribute
				return false;
			}
		}
		else{
			log.error("unknown non-leaf property %s", mainprop.c_str());
			return false;
		}
	}
	else if(input.type == InpType::TIMED_SYNAPTIC){
		// LATER if ever support customization of spike times. Though the number of spikes wouldn't easily change this way, so it's not so useful anyway.
		log.error("timedSynapticInput paths not supported yet");
		return false;
	}
	else{
		// shouldn't happen for now, all types are checked for at this point
		log.error("input source type not supported yet");
		return false;
	}
}
bool Model::ParseLemsQuantityPath_ArtificialCell(const ILogProxy &log, const ArtificialCell &cell,  const std::vector<std::string> &tokens, Simulation::LemsQuantityPath::CellPath &path_cell, Int &tokens_consumed) const {
	if( cell.type == ArtificialCell::SPIKE_SOURCE ){
		path_cell.type = Simulation::LemsQuantityPath::CellPath::INPUT;
		
		const auto &input = input_sources.get( cell.spike_source_seq );
		
		return ParseLemsQuantityPath_InputInstance( log, input, tokens, path_cell.input, tokens_consumed );
	}
	else{
		// work only with LEMSified stuff for now
		if( cell.component.ok() ){
			path_cell.type = Simulation::LemsQuantityPath::CellPath::LEMS;
			return ParseLemsQuantityPathInComponent(log, cell.component, tokens, path_cell.lems_quantity_path, tokens_consumed );
		}
		else{
			log.error("native artificial cell type not supported yet");
			return false;
		}
	}
}

// NOTE: path is a read+write argument. Specifically, segment_seq is read. If set, only the mechanism distributions that include this segment are valid. If unset (ie negative), the caller is responsible for validating that if selected, a mechanism distribution exists for the segments being considered (for example, when multiple segments are considered, is it an error if some segments are outside the distribution mechanism, or are those without the mechanism simply excluded? It depends).
bool Model::ParseLemsQuantityPath_CellProperty(const ILogProxy &log, const CellType &cell_type, const std::vector<std::string> &tokens, Simulation::LemsQuantityPath &path, Int &tokens_consumed ) const {
	
	int segprop_token_id = tokens_consumed;
	
	// now branch according to property type
	if(segprop_token_id >= (Int)tokens.size()){
		log.error("not enough factors for cell or segment property");
		return false;
	}
	const std::string &segprop = tokens[segprop_token_id];
	
	
	if( cell_type.type == CellType::ARTIFICIAL ){
		path.type = Simulation::LemsQuantityPath::CELL;
		return ParseLemsQuantityPath_ArtificialCell(log, cell_type.artificial, tokens, path.cell, tokens_consumed);
	}
	else if(cell_type.type == CellType::PHYSICAL ){
		// use the token
		tokens_consumed++;
		
		const PhysicalCell &cell = cell_type.physical;
		const Morphology &morph = morphologies.get(cell.morphology);
		const BiophysicalProperties &bioph = biophysics.at(cell.biophysicalProperties);
	
		if(segprop == "v"){
			path.type = Simulation::LemsQuantityPath::SEGMENT;
			path.segment.type =  Simulation::LemsQuantityPath::SegmentPath::VOLTAGE;
			//printf("\n voltage yay!\n");
			return true;
		}
		else if(segprop == "caConc"){
			path.type = Simulation::LemsQuantityPath::SEGMENT;
			path.segment.type =  Simulation::LemsQuantityPath::SegmentPath::CALCIUM_INTRA;
			// TODO run a check for existence of Ca2, also for LEMS interface possibly
			return true;
		}
		else if(segprop == "caConc2"){
			path.type = Simulation::LemsQuantityPath::SEGMENT;
			path.segment.type =  Simulation::LemsQuantityPath::SegmentPath::CALCIUM2_INTRA;
			return true;
		}
		else if(segprop == bioph.name || segprop == "biophysicalProperties"){
			
			int biophprop_token_id = segprop_token_id + 1;
			
			if(biophprop_token_id >= (int)tokens.size()){
				log.error("not enough factors for biophysical property");
				return false;
			}
			
			const std::string &bioprop = tokens[biophprop_token_id];
			
			if(bioprop == "membraneProperties"){
				
				// what membrane property is this?
				int biophproprop_token_id = biophprop_token_id + 1;
			
				if(biophproprop_token_id >= (int)tokens.size()){
					log.error("not enough factors for biophysical property");
					return false;
				}
				
				const std::string &bioproprop = tokens[biophproprop_token_id];
				
				// perhaps it is an ion channel distribution?
				Int chandist_seq = -1;
				for( size_t i = 0; i < bioph.membraneProperties.channel_specs.size(); i++ ){
					// NOTE there should be at most one distribution under the same name, TODO discuss and if it's ok refactor to reflect that!
					const auto &dist = bioph.membraneProperties.channel_specs.at(i);
					// if segment_seq is set, select it only if it includes that segment, TODO improve error reporting to "distribution is not present on segment"
					if( dist.name == bioproprop && ((path.segment_seq < 0) || dist.toList(morph).has(path.segment_seq)) ){
						chandist_seq = i;
						break;
					}
					
				}
				if(chandist_seq >= 0){
					// it is
					path.type = Simulation::LemsQuantityPath::CHANNEL;
					path.channel.distribution_seq = chandist_seq;
					const auto &chandist = bioph.membraneProperties.channel_specs[chandist_seq];
					
					int chandistprop_token_id = biophproprop_token_id + 1;
					if(chandistprop_token_id >= (int)tokens.size()){
						log.error("not enough factors for ion distribution property");
						return false;
					}
					const std::string &chandistprop = tokens[chandistprop_token_id];
					
					if(chandistprop == "i"){
						path.channel.type = Simulation::LemsQuantityPath::ChannelPath::I; // all that's needed
						return true;
					}
					else if(chandistprop == "g"){
						path.channel.type = Simulation::LemsQuantityPath::ChannelPath::G; // all that's needed
						return true;
					}
					else if(chandistprop == "iDensity"){
						path.channel.type = Simulation::LemsQuantityPath::ChannelPath::I_DENSITY; // all that's needed
						return true;
					}
					else if(chandistprop == "gDensity"){
						path.channel.type = Simulation::LemsQuantityPath::ChannelPath::G_DENSITY; // all that's needed
						return true;
					}
					else if(chandistprop == ion_channels.getName(chandist.ion_channel)){
						
						// peek inside the channel model
						const auto &channel = ion_channels.get(chandist.ion_channel);
						
						// what channel property is this?
						int ionchanprop_token_id = chandistprop_token_id + 1;
						if(ionchanprop_token_id >= (int)tokens.size()){
							log.error("not enough factors for ion channel property");
							return false;
						}
						const std::string &ionchanprop = tokens[ionchanprop_token_id];
						
						// perhaps it is a gate?
						Int gate_seq = channel.gates.get_id(ionchanprop.c_str());
						if(gate_seq >= 0){
							
							path.channel.gate_seq = gate_seq;
							// what channel gate property is this?
							int gateprop_token_id = ionchanprop_token_id + 1;
							if(gateprop_token_id >= (int)tokens.size()){
								log.error("not enough factors for ion channel gate property");
								return false;
							}
							const std::string &gateprop = tokens[gateprop_token_id];
							
							if(gateprop == "q"){
								path.channel.type = Simulation::LemsQuantityPath::ChannelPath::Q;  // all that's needed
								return true;
							}
							else{
								log.error("unknown ion channel gate property %s", gateprop.c_str());
								return false;
							}
							
						}
						
						//otherwise
						{
							log.error("unknown ion channel property %s", ionchanprop.c_str());
							return false;
						}
						
					}
					else{
						log.error("unknown ion channel distribution property %s", chandistprop.c_str());
						return false;
					}
				}
				// FIXME access all ion pools this way
				
				//otherwise
				{
					log.error("unknown membrane property %s", bioproprop.c_str());
					return false;
				}
				
			}
			else{
				log.error("unknown biophysical property %s", bioprop.c_str());
				return false;
			}
		}
		else if(segprop.find("synapses:") == 0){
			log.error("synapse outputs not supported yet");
			return false;
		}
		else{
			log.error("unknown segment property %s", segprop.c_str());
			return false;
		}
	}
	else{
		log.error("internal error: LEMS quantity path: cell type type %d", cell_type.type);
		return false;
	}
}

bool Model::ParseLemsQuantityPath(const ILogProxy &log, const char *qty_str, const Network &net, Simulation::LemsQuantityPath &path) const {
	auto tokens = string_split(std::string(qty_str), "/");
	const std::string first_identifier = string_split(tokens[0],"[")[0];
	const char *sId = first_identifier.c_str();
	Int tokens_consumed = 0;
	
	// LATER complain if name is ambiguous
	Int group_seq = -1;
	if( (group_seq = net.populations.get_id(sId)) >= 0 ){
		
		if( !ParseLemsSegmentLocator(log, tokens, net, path, tokens_consumed) ) return false;
		// if (cell) segment locator is valid, then it is a property of a cell in a population, and the cell_type is valid.
		const Network::Population &population = net.populations.get(path.population);
		const CellType &cell_type = cell_types.get(population.component_cell);
		
		return ParseLemsQuantityPath_CellProperty(log, cell_type, tokens, path, tokens_consumed );
	}
	else if( (group_seq = net.input_lists.get_id(sId)) >= 0 ){
		// access inputs in an ordered way, it is experimental though
		path.type = Simulation::LemsQuantityPath::INPUT;
		const Network::InputList &list = net.input_lists.get(group_seq);
		const InputSource &input = input_sources.get(list.component);
		// const Int pop_seq = list.population;
		// [](const Network::InputList &group, Int id){return ( 0 <= id && id < (Int)group.input_instances_seq.size() ) ? id : -1 ;}
		path.input.list_seq = group_seq;
		if ( !ParseLemsGroupLocator(log, tokens, "inputList", net.input_lists, [](const Network::InputList &group, Int id){return group.input_instances_seq.getSequential(id);}, path.input.list_seq, path.input.inst_seq, tokens_consumed) ) return false;
		
		// TODO check for weight, delay here maybe?
		return ParseLemsQuantityPath_InputInstance(log, input, tokens, path.input, tokens_consumed);
	}
	else if( (group_seq = net.projections.get_id(sId)) >= 0 ){
		path.type = Simulation::LemsQuantityPath::SYNAPSE;
		const Network::Projection &proj = net.projections.get(group_seq);
		
		if(!ParseLemsGroupLocator(log, tokens, "projection", net.projections, [](const Network::Projection &group, Int id){return group.connections.getSequential(id);}, path.synapse.proj_seq, path.synapse.conn_seq, tokens_consumed)) return false;
		
		// now branch according to property type
		if(!(tokens_consumed+2 <= (Int)tokens.size())){
			log.error("incomplete path for projection element");
			return false;
		}
		
		const Network::Projection::Connection &conn = proj.connections.contents[path.synapse.conn_seq];
		// check for pre, post
		if(!ParseSynapseLocation(log, conn.type, tokens[tokens_consumed].c_str(), path.synapse.location)) return false;
		tokens_consumed++;
		
		Int syn_seq = path.synapse.GetSynSeq(conn);
		const SynapticComponent &syn = synaptic_components.get(syn_seq);
		// TODO check for weight, delay here maybe?
		
		return ParseLemsQuantityPath_SynapticComponent(log, syn, tokens, path.synapse, tokens_consumed);
		// FIXME record a synapse to make sure it works
	}
	// experimental extensions
	else if( (group_seq = net.data_readers.get_id(sId)) >= 0 ){
		path.type = Simulation::LemsQuantityPath::DATAREADER;
		const Network::TimeSeriesReader &reader = net.data_readers.get(group_seq);
		
		path.reader.read_seq = group_seq;
		if(!ParseLemsGroupLocator(log, tokens, "data reader", net.data_readers, [](const auto &reader, Int id){return ( 0 <= id && id < reader.instances ) ? id : -1 ;}, path.reader.read_seq, path.reader.inst_seq, tokens_consumed)) return false;
		// TODO add shortcut without group locator for single instance?
		
		// now branch according to property type
		if(!(tokens_consumed+1 <= (Int)tokens.size())){
			// if the elements have just one property, its name can be skipped
			if(reader.columns.size() == 1){
				path.reader.colu_seq = 0;
				return true;
			}
			else{
				log.error("incomplete path for datareader element");
				return false;
			}
		}
		else{
			const char *sProp = tokens[tokens_consumed].c_str(); tokens_consumed++;
			if(!(tokens_consumed == (Int)tokens.size())){
				log.error("path for datareader element too large");
				return false;
			}
			
			path.reader.colu_seq = reader.columns.get_id(sProp);
			if(path.reader.colu_seq < 0){
				log.error("property %s not found in datareader %s", sProp, net.data_readers.getName(group_seq));
				return false;
			}
			return true;
		}
	}
	// parse them for upcoming refs
	// and get type!
	// FIXME fractionAlong for transmission !
	// TODO flesh this out fully, with special cases, ions, and such.
	// TODO make a reverse toString() for lemsquantitypaths.
	// path locations: parse, gettype, print, to tabloc.
	else{
		log.error("unknown top level identifier \"%s\" for quantity path", sId);
		return false;
	}
	
	// TODO access projections and input_lists.
}

bool Model::ParseLemsEventPathInComponent( const ILogProxy &log, const ComponentInstance &instance, const std::vector<std::string> &tokens, Simulation::LemsInstanceEventPath &liep, Int &tokens_consumed) const {
	
	if(tokens_consumed >= (Int)tokens.size()){
		log.error("not enough factors for LEMS component event path");
		return false;
	}
	
	// nothing more to extract
	const char *sEventPort = tokens[tokens_consumed].c_str();
	tokens_consumed++;
	
	if( tokens_consumed < (int)tokens.size() ){
		log.error( "LEMS child component event outputs not yet supported" );
		return false;
	}
	
	// get comp.type
	Int comptype_seq = instance.id_seq;
	if( !component_types.has(comptype_seq) ){
		log.error("internal error: LEMS event path missing component type %d", (int) instance.id_seq);
		return false;
	}
	
	// which takes preference in namespace anyway?
	const auto &comptype = component_types.get(comptype_seq);
	
	Int &event_port_seq = liep.event_port_seq = -1;
	Int event_in_seq  = comptype.event_inputs .get_id(sEventPort);
	Int event_out_seq = comptype.event_outputs.get_id(sEventPort);
	if( event_in_seq < 0 && event_out_seq < 0 ){
		log.error("%s is not a defined event port in component type %s", sEventPort, component_types.getName(comptype_seq));
		return false;
	}
	else if( event_in_seq >= 0 && event_out_seq >= 0 ){
		log.error("%s is both an input and output port", sEventPort);
		return false;
	}
	else if( event_in_seq >= 0 ){
		liep.type = Simulation::LemsInstanceEventPath::IN;
		event_port_seq = event_in_seq;
		
		// XXX this is just a limitation in Eden right now
		// but idk if input event ports are logged in practice
		log.error("%s is an input event port; which is not yet supported in EDEN", sEventPort);
		return false;
	}
	else if( event_out_seq >= 0 ){
		liep.type = Simulation::LemsInstanceEventPath::OUT;
		event_port_seq = event_out_seq;
		
		return true;
	}
	else{
		// all cases are handled but the compiler whines for some reason
		assert(false);
		return false;
	}
}
bool Model::ParseLemsEventPath_InputInstance( const ILogProxy &log, const InputSource &input, const std::vector<std::string> &tokens, Simulation::InputInstanceEventPath &path, Int &tokens_consumed ) const {
	typedef InputSource::Type InpType;
	typedef Simulation::InputInstanceEventPath InpPath;
	int mainprop_token_id = tokens_consumed;
	
	// now branch according to property type
	if(mainprop_token_id >= (Int)tokens.size()){
		log.error("not enough factors for input event path");
		return false;
	}
	const std::string &mainprop = tokens[mainprop_token_id];
	
	// const bool is_native = ( input.type == InpType::PULSE );
	const bool has_synapse = (
		input.type == InpType::TIMED_SYNAPTIC
		|| input.type == InpType::POISSON_SYNAPSE
		|| input.type == InpType::POISSON_SYNAPSE_TRANSIENT
	);
	
	// synapse containing stuff are more than just lemsified
	if( has_synapse && mainprop == "synapse" ){
		tokens_consumed++;
		log.error("event paths in synpase children of inputs not supported yet");
		return false;
	}
	
	// work only with LEMSified stuff for now
	if( input.component.ok() ){
		path.type = InpPath::LEMS;
		return ParseLemsEventPathInComponent(log, input.component, tokens, path.lems_event_path, tokens_consumed );
	}
	// not lemsified is propably not supported
	// but let's check if there was any spike port even exists in native components (not yet)
	
	// NB: not looking for other children yet ...
	if(!( input.HasSpikeOut(component_types))){ // || input.HasSpikeIn(component_types)
		log.error("input source type has no event ports");
		return false;
	}
	
	// otherwise it's a native, or complicated, type.
	path.type = InpPath::NATIVE;
	
	tokens_consumed++; // note: this token will be used in the following, unshift or refactor if necessary later
	
	auto &native_entry = path.native_entry;
	typedef InpPath::NativeEntry Nat;
	
	// is it a terminal property? (to avoid repeated checks in the following)
	if(tokens_consumed == (Int)tokens.size()){
		if(     mainprop == "spike"){native_entry = Nat::SPIKE; return true; }
		else{
			log.error("unknown leaf property %s", mainprop.c_str()); // LATER log in token possibly in attribute
			return false;
		}
	}
	else{
		log.error("unknown non-leaf property %s", mainprop.c_str());
		return false;
	}
	
	if( input.type == InpType::SPIKE_LIST ){
		if(( input.type == InpType::SPIKE_LIST && !( native_entry == Nat::SPIKE ) )
		){
			log.error("property %s does not match input source type", mainprop.c_str()); return false;
		}
		return true; // all checks passed
	}
	else{
		// shouldn't happen for now, all types are checked for at this point
		log.error("input source type not supported yet");
		return false;
	}
}
bool Model::ParseLemsEventPath_ArtificialCell(const ILogProxy &log, const ArtificialCell &cell,  const std::vector<std::string> &tokens, Simulation::LemsEventPath::Cell &path_cell, Int &tokens_consumed) const {	
	if( cell.type == ArtificialCell::SPIKE_SOURCE ){
		path_cell.type = Simulation::LemsEventPath::Cell::INPUT;
		
		const auto &input = input_sources.get( cell.spike_source_seq );
		
		return ParseLemsEventPath_InputInstance( log, input, tokens, path_cell.input, tokens_consumed );
	}
	else{
		// work only with LEMSified stuff for now
		if( cell.component.ok() ){
			path_cell.type = Simulation::LemsEventPath::Cell::LEMS;
			return ParseLemsEventPathInComponent(log, cell.component, tokens, path_cell.lems_event_path, tokens_consumed );
		}
		else{
			log.error("native artificial cell type not supported yet");
			return false;
		}
	}
}
// NOTE: path is a read+write argument. Specifically, segment_seq is read. If set, only the mechanism distributions that include this segment are valid. If unset (ie negative), the caller is responsible for validating that if selected, a mechanism distribution exists for the segments being considered (for example, when multiple segments are considered, is it an error if some segments are outside the distribution mechanism, or are those without the mechanism simply excluded? It depends).
bool Model::ParseLemsEventPath_CellProperty(const ILogProxy &log, const CellType &cell_type, const std::vector<std::string> &tokens, Simulation::LemsEventPath &path, Int &tokens_consumed ) const {

	int segprop_token_id = tokens_consumed;
	
	// now branch according to property type
	if(segprop_token_id >= (Int)tokens.size()){
		log.error("not enough factors for cell or segment event path");
		return false;
	}
	
	if( cell_type.type == CellType::ARTIFICIAL ){
		path.type = Simulation::LemsEventPath::CELL;
		return ParseLemsEventPath_ArtificialCell(log, cell_type.artificial, tokens, path.cell, tokens_consumed);
	}
	else if(cell_type.type == CellType::PHYSICAL ){
		const std::string &segprop = tokens[segprop_token_id];
		segprop_token_id++;
		if( segprop_token_id == (Int)tokens.size() ){
			// right on the segment
			// TODO check if Vthreshold is defined for that segment
			
			std::string port = segprop;
			if(port == "spike"){
				path.type = Simulation::LemsEventPath::SEGMENT;
				path.segment.type =  Simulation::LemsEventPath::Segment::SPIKE;
				//printf("\n spikes yay!\n");
				return true;
			}
			else{
				log.error("unknown eventPort %s", port.c_str());
				return false;
			}
		}
		else{
			printf("%d, %d, %s!%s|%s|%s\n", segprop_token_id, (int)tokens.size(), segprop.c_str(), tokens[0].c_str(), tokens[1].c_str(), tokens[2].c_str());
			// will keep on parsing the path LATER, if ever
			log.error("spiking subcomponents of neuron segment not supported yet");
			return false;
		}
	}
	else{
		log.error("internal error: LEMS event path: cell type type %d", cell_type.type);
		return false;
	}
}
// TODO merge into a common skeleton with quantity path parsing
bool Model::ParseLemsEventPath(const ILogProxy &log, const char *sPath, const char *port_name_or_null_or_empty, const Network &net, Simulation::LemsEventPath &path) const {
	auto tokens = string_split(std::string(sPath), "/");
	if( port_name_or_null_or_empty && *port_name_or_null_or_empty ) tokens.push_back(port_name_or_null_or_empty);
	const std::string first_identifier = string_split(tokens[0],"[")[0];
	const char *sId = first_identifier.c_str();
	Int tokens_consumed = 0;
	
	// LATER complain if name is ambiguous
	Int group_seq = -1;
	if( (group_seq = net.populations.get_id(sId)) >= 0 ){
	
		if ( !ParseLemsSegmentLocator(log, tokens, net, path, tokens_consumed) ) return false;
		// if (cell) segment locator is valid, then it is a property of a cell in a population, and the cell_type is valid.
		const Network::Population &population = net.populations.get(path.population);
		const CellType &cell_type = cell_types.get(population.component_cell);
		
		return ParseLemsEventPath_CellProperty(log, cell_type, tokens, path, tokens_consumed );
	}
	// experimental extensions
	else if( (group_seq = net.event_readers.get_id(sId)) >= 0 ){
		path.type = Simulation::LemsEventPath::EVENTREADER;
		const Network::EventSetReader &reader = net.event_readers.get(group_seq);
		path.reader.read_seq = group_seq;
		if(!ParseLemsGroupLocator(log, tokens, "event reader", net.event_readers, SeqToSeq(&Network::EventSetReader::instances), path.reader.read_seq, path.reader.inst_seq, tokens_consumed)) return false;
		// TODO add shortcut without group locator for single instance?
		
		// now branch according to property type
		if(!(tokens_consumed+1 <= (Int)tokens.size())){
			// if the elements have just one property, its name can be skipped
			if(reader.ports.size() == 1){
				path.reader.port_seq = 0;
				return true;
			}
			else{
				log.error("incomplete path for eventreader element");
				return false;
			}
		}
		else{
			const char *sProp = tokens[tokens_consumed].c_str(); tokens_consumed++;
			if(!(tokens_consumed == (Int)tokens.size())){
				log.error("path for eventreader element too large");
				return false;
			}
			
			path.reader.port_seq = reader.ports.get_id(sProp);
			if(path.reader.port_seq < 0){
				log.error("property %s not found in eventreader %s", sProp, net.event_readers.getName(group_seq));
				return false;
			}
			return true;
		}	
	}
	else{
		log.error("unknown top level identifier \"%s\" for event path", sId);
		return false;
	}
}

//------------------> Debug printing

void AcrossSegOrSegGroup::debug_print(const Morphology &morph) const{
	if(type == SEGMENT) printf("%ld", seqid);
	else if(type == GROUP) printf("%s", morph.Stringify(morph.segment_groups[seqid]).c_str() );
	else printf("what???");
	printf("\n");
}

std::string SpeciesAcrossSegOrSegGroup::Stringify(const CollectionWithNames<IonSpecies> &ion_species) const {
	char buf[300];
	sprintf(buf," initial int: %g ext: %g %s", initialConcentration, initialExtConcentration, Scales<Concentration>::native.name);
	return "Ion species #" +std::to_string(species)+buf;
}

void BiophysicalProperties::debug_print(const Morphology &morph, const CollectionWithNames<IonSpecies> &ion_species) const {
	//debug output
	printf("Biophysics contents\n");
	
	printf("Intracellular properties:\n");
	for(auto spec : intracellularProperties.resistivity_specs){
		printf("Resistivity: %g %s for ", spec.value, Scales<Resistivity>::native.name);
		spec.debug_print(morph);
	}
	
	for(auto spec : intracellularProperties.ion_species_specs){
		printf("%s for ", spec.Stringify(ion_species).c_str() );
		spec.debug_print(morph);
	}
	//TODO extract name to display species
	printf("\n");
	
	
	if(!extracellularProperties.ion_species_specs.empty()){
		printf("Extracellular properties:\n");
		for(auto spec : extracellularProperties.ion_species_specs){
			printf("%s for ", spec.Stringify(ion_species).c_str() );
			spec.debug_print(morph);
		}
	}
	printf("\n");
	
	printf("Membrane properties:\n");
	for(auto spec : membraneProperties.capacitance_specs){
		printf("Specific capacitance: %g %s for ", spec.value, Scales<SpecificCapacitance>::native.name);
		spec.debug_print(morph);
	}
	for(auto spec : membraneProperties.initvolt_specs){
		printf("Initial voltage: %g %s for ", spec.value, Scales<Voltage>::native.name);
		spec.debug_print(morph);
	}
	for(auto spec : membraneProperties.threshold_specs){
		printf("Spike threshold: %g %s for ", spec.value, Scales<Voltage>::native.name);
		spec.debug_print(morph);
	}
	
}

void ComponentType::debug_print( const DimensionSet &dimensions ) const {
	
	for(auto keyval : constants.names){
		const auto &thing = constants.get(keyval.second);
		printf( "Constant %s: %s %f", keyval.first, dimensions.Stringify(thing.dimension).c_str(), thing.value );
		if(thing.dimension != Dimension::Unity()) printf( " (%s)",  dimensions.GetNative(thing.dimension).name.c_str() );
		printf("\n");
	}
	
	for(auto keyval : properties.names){
		const auto &thing = properties.get(keyval.second);
		printf( "Parameter %s: %s %f", keyval.first, dimensions.Stringify(thing.dimension).c_str(), thing.value );
		if(thing.dimension != Dimension::Unity()) printf( " (%s)",  dimensions.GetNative(thing.dimension).name.c_str() );
		printf("\n");
	}
	
	for(auto keyval : exposures.names){
		printf("Exposure %s: %s\n", keyval.first, dimensions.Stringify( getExposureDimension(keyval.first) ).c_str());
	}
	
	for(auto keyval : requirements.names){
		const auto &thing = requirements.get(keyval.second);
		printf("Requirement %s: %s\n", keyval.first, dimensions.Stringify(thing.dimension).c_str());
	}
	
	for(auto keyval : event_inputs.names){
		const auto &thing = event_inputs.get(keyval.second);
		printf("Event input %s\n", keyval.first);
		(void) thing;		
	}
	for(auto keyval : event_outputs.names){
		const auto &thing = event_outputs.get(keyval.second);
		printf("Event output %s\n", keyval.first);
		(void) thing;
	}
	printf("\n");
};

//------------------> NeuroML Parsing

const char * RequiredNmlId(const ImportLogger &log, const pugi::xml_node &from_node){
	const char *ret_id = from_node.attribute("id").value(); //if missing, null handle returns empty string
	if(!*ret_id){
		log.error(from_node, "element lacks required NML ID");
		return NULL;
	}
	return ret_id;
}
// a variant for the common CollectionWithNames
template<typename T>
const char * RequiredNmlId(const ImportLogger &log, const pugi::xml_node &from_node, const CollectionWithNames<T> &set, const char *parent_name = NULL){
	const char *name = RequiredNmlId(log, from_node);
	if(!name) return NULL;
	if(set.has(name)){
		if(name) log.error(from_node, "%s %s already defined in %s", from_node.name(), name, parent_name);
		else log.error(from_node, "%s %s already defined", from_node.name(), name);
		return NULL;
	}
	return name;
}
// just in case it may differ more than the name, and to discern between NeuroML and pure LEMS
const char * RequiredLemsName(const ImportLogger &log, const pugi::xml_node &from_node){
	const char *ret_id = from_node.attribute("name").value(); //if missing, null handle returns empty string
	if(!*ret_id){
		log.error(from_node, "element lacks required LEMS name");
		return NULL;
	}
	return ret_id;
}
// TODO find and replace with if(!*
const char * RequiredAttribute(const ImportLogger &log, const pugi::xml_node &from_node, const char *attr_name){
	const char *ret_attr = from_node.attribute(attr_name).value(); //if missing, null handle returns empty string
	if(!*ret_attr){
		log.error(from_node, "must have \"%s\" attribute", attr_name);
		return NULL;
	}
	return ret_attr;
}

bool RequiredPositiveAttribute(const ImportLogger &log, const pugi::xml_node &from_node, const char *attr_name, Int &value){
	const char *sInt = RequiredAttribute( log, from_node, attr_name);
	if(!sInt) return false;
	if(!( StrToL(sInt, value) && value > 0 )){
		log.error(from_node, " \"%s\" must be a positive integer, not %s", attr_name, sInt);
		return false;
	}
	return true;
}
bool RequiredNonnegativeAttribute(const ImportLogger &log, const pugi::xml_node &from_node, const char *attr_name, Int &value){
	const char *sInt = RequiredAttribute( log, from_node, attr_name);
	if(!sInt) return false;
	if(!( StrToL(sInt, value) && value >= 0 )){
		log.error(from_node, " \"%s\" must be a non-negative integer, not %s", attr_name, sInt);
		return false;
	}
	return true;
}
bool OptionalPositiveAttribute(const ImportLogger &log, const pugi::xml_node &from_node, const char *attr_name, Int &value){
	if(*from_node.attribute(attr_name).value()) return RequiredPositiveAttribute(log, from_node, attr_name, value);
	else return true;
}
bool OptionalNonnegativeAttribute(const ImportLogger &log, const pugi::xml_node &from_node, const char *attr_name, Int &value){
	if(*from_node.attribute(attr_name).value()) return RequiredNonnegativeAttribute(log, from_node, attr_name, value);
	else return true;
}

const char * RequiredComponentType(const ImportLogger &log, const pugi::xml_node &from_node){
	const char *sType = from_node.name();
	if( strcmp(sType, "Component") == 0 ){
		sType = from_node.attribute("type").value();
		if(!*sType){
			log.error(from_node, "<Component> must have a \"type\" attribute");
			return NULL;
		}
		// printf("actually a %s\n", sType);
	}
	assert( *sType ); // assume tag name is not empty
	return sType;
}
// TODO refactor this better, consider RequiredComponentType instead
const char * TagNameOrType(const pugi::xml_node &from_node){
	auto ret = from_node.attribute("type").value();
	if( *ret ) return ret;
	else return from_node.name();
}

// TODO if it's a file url, much like the current behaviour of OutputFile, it's relative to the cwd of the program, not the file being parsed! TODO specify a way to selecthow to behave, or at least document --  filename
// here's a trick: use "no uri scheme" to mean "relative to xml file" and "file://" to mean "relative to working directory", it is also how href's in webpages work, TODO document it.
bool GetSchemeOrRelpathFromUrl(const ImportLogger &log, const pugi::xml_node &from_node, const char *attrname,  std::string &url ){
	
	const char *sHref = RequiredAttribute( log, from_node, attrname);
	if(!sHref) return false;
	
	url = sHref;
	
	std::string scheme, auth_path;
	if(!GetUrlScheme(url, scheme, auth_path)){
		const char *loading_from_file = log.GetFilenameFromElement(from_node);
		url = GetRelativeFilePath((loading_from_file ? loading_from_file : "."), url);
	}
	// else keep the url as is
	
	return true;
}

bool ParseMorphology(const ImportLogger &log, const pugi::xml_node &eMorph, Morphology &morph ){
	
	//contains segments and segment groups
	
	// keep references to tags of segments, for error reporting
	// CollectionWithIds<const pugi::xml_node *> tags_per_segment;
	
	//Note that strict segment id ordering and parent < child ensure loops cannot be formed
	for (auto eMorphEl: eMorph.children()){
		// printf("elm  %s...\n", eMorphEl.name());
		// perhaps annotation stuff?
		// NOTE wherever there is Standalone, notes, annotation and property tage may exist. Handle properties as they appear, for they must be documented to be used.
		if(
			strcmp( eMorphEl.name(), "notes" ) == 0
			|| strcmp( eMorphEl.name(), "annotation" ) == 0
		){
			continue;
		}
		else if(strcmp(eMorphEl.name(), "segment") == 0){
			const auto &eSeg = eMorphEl;
			Morphology::Segment seg;
			
			//get required id and check if it is in typical id order
			Int nml_seg_id = -1;
			if( !StrToL(eSeg.attribute("id").value(), nml_seg_id) ){
				if(!eSeg.attribute("id")) log.error(eSeg, "segment tag lacks id attribute");
				else if(errno == ERANGE) log.error(eSeg, "segment id %s out of range", eSeg.attribute("id").value());
				else log.error(eSeg, "segment id %s not a number", eSeg.attribute("id").value());
				return false;
			}
			if(nml_seg_id < 0){
				log.error(eSeg, "segment id is negative");
				return false;
			}
			if(morph.segments.hasId(nml_seg_id)){
				log.error(eSeg, "segment %s already defined", eSeg.attribute("id").value());
				return false;
			}
			// printf("elm  %ld...\n", nml_seg_id);
			//If parent element exists, get id and optionally fractionAlong
			Int nml_seg_parent;
			seg.parent = -1;
			seg.fractionAlong = 1;
			auto eSegParent = eSeg.child("parent");
			if(eSegParent){
				if( !StrToL(eSegParent.attribute("segment").value(), nml_seg_parent) ){
					log.error(eSegParent, "segment %ld has an invalid parent id %s", nml_seg_id, eSegParent.attribute("segment").value());
					return false;
				}
				//should be negative(ie null) or an already existing segment
				if(nml_seg_parent < 0){
					seg.parent = -1;
				}
				else{
					Int internal_seg_parent = morph.lookupSegId(nml_seg_parent);
					if( internal_seg_parent < 0 ){
						log.error(eSegParent, "segment %ld is child of not previously defined segment %ld", nml_seg_id, nml_seg_parent);
						return false;
					}
					seg.parent = internal_seg_parent;
				}
				
				
				if(eSegParent.attribute("fractionAlong")){
					if(!(
						StrToF(eSegParent.attribute("fractionAlong").value(), seg.fractionAlong)
						&& 0 <= seg.fractionAlong && seg.fractionAlong <= 1
					)){
						log.error(eSegParent, "segment %ld has an invalid fractionAlong %s", nml_seg_id, eSegParent.attribute("fractionAlong").value());
						return false;
					}
				}
			}
			// printf("distal  %ld...\n", nml_seg_id);
			
			//get distal part of segment
			auto eSegDistal = eSeg.child("distal");
			if(!eSegDistal){
				log.error(eSeg, "segment %ld has no distal point", nml_seg_id);
				return false;
			}
			
			if(!(
				   StrToF(eSegDistal.attribute("x").value(), seg.distal.x)
				&& StrToF(eSegDistal.attribute("y").value(), seg.distal.y)
				&& StrToF(eSegDistal.attribute("z").value(), seg.distal.z)
				&& StrToF(eSegDistal.attribute("diameter").value(), seg.distal.d)
				&& seg.distal.d >= 0
			)){
				log.error(eSegDistal,
					"segment %ld has an invalid distal point (%s, %s, %s), %s (x,y,z),diameter combination", nml_seg_id,
					eSegDistal.attribute("x").value(),
					eSegDistal.attribute("y").value(),
					eSegDistal.attribute("z").value(),
					eSegDistal.attribute("diameter").value()
				);
				return false;
			}
			// distal diameter might be allowed to be zero, this is checked along with proximal diameter
			
			// printf("proximal  %ld...\n", nml_seg_id);
			//get proximal part of segment, or interpolate from parent
			auto eSegProximal = eSeg.child("proximal");
			if(!eSegProximal){
				//interpolate precise attachment point to the parent segment
				if(seg.parent < 0){	
					log.error(eSeg, "root segment %ld has no proximal point defined", nml_seg_id);
					return false;
				}
				// printf("parent  %ld...\n", seg.parent);
				const auto &parent_seg = morph.segments.atSeq(seg.parent);
				
				Real r = seg.fractionAlong; const auto &p = parent_seg.proximal; const auto &d = parent_seg.distal;
				seg.proximal.x = Lerp(p.x, d.x, r);
				seg.proximal.y = Lerp(p.y, d.y, r);
				seg.proximal.z = Lerp(p.z, d.z, r);
				seg.proximal.d = Lerp(p.d, d.d, r);
				
			}else{
				//explicitly read it
				if(!(
					   StrToF(eSegProximal.attribute("x").value(), seg.proximal.x)
					&& StrToF(eSegProximal.attribute("y").value(), seg.proximal.y)
					&& StrToF(eSegProximal.attribute("z").value(), seg.proximal.z)
					&& StrToF(eSegProximal.attribute("diameter").value(), seg.proximal.d)
					&& seg.proximal.d >= 0
				)){
					log.error(eSegProximal,
						"segment %ld has an invalid proximal point (%s, %s, %s), %s (x,y,z),diameter combination", nml_seg_id,
						eSegProximal.attribute("x").value(),
						eSegProximal.attribute("y").value(),
						eSegProximal.attribute("z").value(),
						eSegProximal.attribute("diameter").value()
					);
					return false;
				}
				// proximal diameter might be allowed to be zero, this is checked along with distal diameter and parent relationship
			}
			
			// check that no intermediate segment is connected at/with a point of zero diameter
			if( seg.proximal.d <= 0 && seg.distal.d <= 0 ){
				log.error(eSeg, "segment %ld cannot have proximal and distal edge diameters both zero", nml_seg_id);
				return false;
			}
			// if the segment is attached to a parent, proximal.d must be non-zero
			if( seg.parent != 1){
				// proximal.d = 0 could be caused either explicitly, or by conneting to a parent that tapers to zero on its fractionAlong = 1
				if( seg.proximal.d <= 0 ){
					log.error(eSeg, "segment %ld cannot have a zero-diameter contact point with parent", nml_seg_id);
				}
			}
			// whether distal.d = 0 is allowed will be checked when a child connects to it. FIXME make sure Eden avoids integrating the taper to zero if fractionAlong != 1
			
			// also keep a reference to this element for further validation, if needed
			// tags_per_segment.add(&eSeg, nml_seg_id);
			
			//add new segment to Morphology array, yay!
			morph.addSegment(seg, nml_seg_id);
		}
		else if(strcmp(eMorphEl.name(), "segmentGroup") == 0){
			const auto &eGroup = eMorphEl;
			
			Morphology::SegmentGroup group;
			
			//could also dump all segments and stabilize in the end, for more efficiency
			
			//must have a string ID, so it can be referred to
			const char * group_name = eGroup.attribute("id").value();
			if(!( *group_name )){
				//name is empty or missing
				log.error(eGroup, "group has no name");
				return false;
			}
			if(morph.segment_groups_by_name.count(group_name)){
				log.error(eGroup, "group %s already defined", group_name);
				return false;
			}
			
			const char * sNeuroLexId = eGroup.attribute("neuroLexId").value();
			if(strcmp(sNeuroLexId, "sao864921383") == 0){
				// This mystery tag affects how a section of neurite is discretised !
				// The segments in the group must form a continuous non-branching section.
				// The segments are merged in one simulated compartment by default.
				// If the <property tag="numberInternalDivisions" value="(number here)"/> tag is present in the group,
				// the section is discretised along equal central path lengths, as NEURON does for pt3d based descriptions of sections.
				// (Refer to the NEURON book chapter 5  for more details, and figure 5.6. for an example.)
				// Note that this could be a suboptimal discretisation method if the neurite diameter, or parameters of active mechanisms, vary wildly. Take care when building models!
				group.is_cable = true;
				
				// validate its being a cable after all child elements are parsed
			}
			// there are more commonly named neuroLexId's like "soma" "axon" "dendrites" but they don't affect the model
			
			for (auto eGroupEl: eGroup.children()){
				if(strcmp(eGroupEl.name(), "member") == 0){
					
					//add one segment to group
					Int nml_memseg = -1;
					if( !StrToL(eGroupEl.attribute("segment").value(), nml_memseg) ){
						log.error(eGroupEl, "group %s has an invalid member id %s", group_name, eGroupEl.attribute("segment").value());
						return false;
					}
					Int internal_memseg = morph.lookupSegId(nml_memseg);
					if( internal_memseg < 0 ){
						log.error(eGroupEl, "group %s has a not previously defined segment %ld", group_name, nml_memseg);
						return false;
					}
					
					group.addd(internal_memseg);
					
				}
				else if(strcmp(eGroupEl.name(), "include") == 0){
					//add one group to group
					const char *memgroup = eGroupEl.attribute("segmentGroup").value();
					if(!( *memgroup )){
						//name is empty or missing
						log.error(eGroupEl, "included segment group has no name");
						return false;
					}
					if(!morph.segment_groups_by_name.count(memgroup)){
						log.error(eGroup, "included segment group %s not already defined", memgroup);
						return false;
					}
					Int memid = morph.segment_groups_by_name[memgroup];
					// printf("add seggroup %s to %s \n", memgroup, group_name);
					const auto &addgroup = morph.segment_groups[memid];
					group.add( addgroup );
				}
				else if(strcmp(eGroupEl.name(), "path") == 0){
					// add path from segment to segment
					// (may include soma i guesss? what if later non-tree topologies are supported?)
					log.error(eGroupEl, "group %s has path, not supported yet", group_name);
					return false;
				}
				else if(strcmp(eGroupEl.name(), "subTree") == 0){
					//add one group to group
					log.error(eGroupEl, "group %s has subTree, not supported yet", group_name);
					return false;
				}
				else if(strcmp(eGroupEl.name(), "inhomogeneousParameter") == 0){
					//add one metric to group
					// actually, add it to entire morphology because why not
					//see SubsetDomainIterator hoc in NEURON Exporter
					//parameter is hardcoded to be g_max, stuff is still under construction
					const auto &eInhoParm = eGroupEl;
					
					Morphology::InhomogeneousParameter inhoparm;
					
					auto name = eInhoParm.attribute("id").value();
					if( !*name ){
						log.error(eInhoParm, "id attribute missing");
						return false;
					}
					if( group.inhomogeneous_parameters.has(name) ){
						log.error(eInhoParm, "id %s already defined in segment group", name);
						return false;
					}
					
					auto variable = eInhoParm.attribute("variable").value();
					if( !*variable ){
						log.error(eInhoParm, "variable attribute missing");
						return false;
					}
					// could also verify that it doesn't have the same 'variable name' as another variable, that is best handled in the specific case that might bring two variables together in the future
					inhoparm.variable = variable;
					
					auto metric = eInhoParm.attribute("metric").value();
					if( !*metric ){
						log.error(eInhoParm, "metric attribute missing");
						return false;
					}
					if( stricmp(metric, "Path Length from root") != 0 ){
						log.error(eInhoParm, "unknown metric %s", metric);
						return false;
					}
					inhoparm.metric = Morphology::InhomogeneousParameter::PATH_LENGTH_FROM_ROOT;
					
					// note: pugixml returns "" if the child, attr etc. does not exist (thus a null handle is used)
					// Break the access down for alterantive parsers later
					auto eProximal = eInhoParm.child("proximal");
					auto translationStart = eProximal.attribute("translationStart").value();
					
					if(!*translationStart){
						inhoparm.subtract_the_minimum = false;
					}
					else{
						Real val;
						if(!(
							StrToF(translationStart, val)
							&& val == 0
						)){
							log.error(eProximal, "invalid %s value \"%s\"", "translationStart", translationStart);
							return false;
						}
						
						inhoparm.subtract_the_minimum = true;
					}
					auto eDistal = eInhoParm.child("distal");
					auto normalizationEnd = eDistal.attribute("normalizationEnd").value();
					
					if(!*normalizationEnd){
						inhoparm.divide_by_maximum = false;
					}
					else{
						Real val;
						if(!(
							StrToF(normalizationEnd, val)
							&& val == 1
						)){
							log.error(eProximal, "invalid %s value \"%s\"", "normalizationEnd", normalizationEnd);
							return false;
						}
						
						inhoparm.divide_by_maximum = true;
					}
					
					// add inhomogeneous parameter to segment group, yay!
					group.inhomogeneous_parameters.add(inhoparm, name);
				}
				else if(strcmp(eGroupEl.name(), "property") == 0){
					auto sTag = eGroupEl.attribute("tag").value();
					auto sValue = eGroupEl.attribute("value").value();
					
					// NeuroML v2 would prefer each segment to map to a single abstraction entity,
					// but v1 and NEURON models had to be imported too, so that's a compatibility/reproducibility feature
					// but yet there is a real need to represent a section with many 3d points and few compartments
					if( strcmp(sTag, "numberInternalDivisions") == 0 ){
						if( !group.is_cable ){
							log.error(eGroupEl, "group %s has attribute %s but is not a cable (add neuroLexId=\"sao864921383\" attribute to <segmentGroup> to make the group a cable)", group_name, sTag);
							return false;
						}
						
						if(!(
							StrToL(sValue, group.cable_nseg)
							&& group.cable_nseg > 0
						)){
							log.error(eGroupEl, "group %s has an invalid %s parameter : %s", group_name, sTag, sValue);
							return false;
						}
						// cable_nseg has a positive integer value, that is enough validation
					}
					else{
						log.error(eGroupEl, "unknown property %s", sTag);
						return false;
					}
				}
				else{
					//fall through
				}
			}
			
			//compact the group after all these unions
			group.list.Compact();
			if(group.is_cable){
				// validate its being a cable on the spot
				// due to the order of appearance of segment, the ultimate descendant must be the last one
				// and its sequence of ancestors must be the group id's in reverse order
				std::vector<Int> segs = group.list.toArray();
				for(int i = segs.size() - 1; i > 0; i--){
					const auto &seg = morph.segments.atSeq(segs[i]);
					// check that the parent is correct
					if(seg.parent != segs[i-1]){
						log.error(eGroup, "group is not a proper unbranched cable: segment %ld 's parent is %ld when it should be %ld",
							morph.segments.getId(segs[i]), morph.segments.getId(seg.parent), morph.segments.getId(segs[i-1]));
						return false;
					}
					// check that fractionAlong = 1 (attached to end of parent) for a properly unbranched cable
					if(seg.parent >= 0 && seg.fractionAlong != 1){
						log.error(eGroup, "group is not a proper unbranched cable: segment %ld is not attached to end of parent but on fractionAlong = %.17g instead",
							morph.segments.getId(segs[i]), seg.fractionAlong);
						return false;
					}
				}
				
				// validate its being non-overlapping after all groups are parsed
			}
			
			//add group to Morphology, yay!
			morph.add(group, group_name);
		}
		else{
			log.error( eMorphEl, "unknown morphology tag %s", eMorphEl.name() );
			return false;
		}
		
	} //end morph children
	
	//if there is no "all" group, add one so it will be available for future reference
	if(!morph.segment_groups_by_name.count("all")){
		Morphology::SegmentGroup group;
		group.list = morph.getFullList();
		morph.add(group, "all");
	}
	
	//ensure there is at least one segment
	if(morph.segments.contents.empty()){
		log.error(eMorph, "there should be at least one segment");
		return false;
	}
	// verify that the cables do not overlap
	bool morph_uses_cables = false;
	std::vector<Int> cable_per_segment(morph.segments.size(), -1); //group_seq of cable for each segment
	
	for( Int group_seq = 0; group_seq < (Int) morph.segment_groups.size(); group_seq++){
		const auto &group = morph.segment_groups[group_seq];
		
		if( !group.is_cable ) continue;
		morph_uses_cables = true;
		
		std::vector<Int> segs = group.list.toArray();
		for( Int seg_seq : segs ){
			
			// perhaps the semgnet is already included in another group?
			Int previous_group_seq = cable_per_segment[seg_seq];
			if(!(previous_group_seq < 0)){
				// cables should not overlap!!
				
				// FIXME show names!
				// log.error(eMorph, "segment groups %s and %s should be non-overlapping, yet they overlap on segment %ld", , , morph.lookupNmlId(seg_seq));
				log.error(eMorph, "segment groups %ld and %ld should be non-overlapping, yet they overlap on segment %ld", group_seq, previous_group_seq, morph.lookupNmlId(seg_seq));
				// maybe show offending group tags LATER?
				return false;
			}
			cable_per_segment[seg_seq] = group_seq;
		}
		
	}
	
	// if cables exist and yet they do not form the whole cell, warn
	if(morph_uses_cables){
		IdListRle segs_not_covered_by_cables;
		for( Int seg_seq = 0; seg_seq < (Int) cable_per_segment.size(); seg_seq++ ){
			if(cable_per_segment[seg_seq] < 0) segs_not_covered_by_cables.Addd(seg_seq);
		}
		segs_not_covered_by_cables.Compact();
		if( segs_not_covered_by_cables.Count() > 0 ){
			log.warning(eMorph, "the morphology's non-overlapping cables do not form the whole morphology. Segments not included in cables: ", morph.Stringify_SegSeq_List(segs_not_covered_by_cables).c_str());
		}
	}
	
	if(log.debug_printing) morph.debug_print();
	
	return true;
}

bool ParseAcrossSegGroup( const ImportLogger &log, const char *sGroupName, const pugi::xml_node &eAppliedOn, const Morphology &morph, AcrossSegOrSegGroup &applied_on ){
	if(!sGroupName || !*sGroupName){
		//check for "all" segmentGroup
		//could possibly also create an actual "all" group, but who knows how "all" is actually used
		sGroupName = "all";
	}
	//check if it exists in morphology
	if(!morph.segment_groups_by_name.count(sGroupName)){
		log.error(eAppliedOn, "group %s does not exist in associated Morphology", sGroupName);
		return false;
	}
	
	applied_on.Group(morph.segment_groups_by_name.at(sGroupName));
	return true;
}
bool ParseAcrossSegOrSegGroup(const ImportLogger &log, const pugi::xml_node &eAppliedOn, const Morphology &morph, AcrossSegOrSegGroup &applied_on){
	//Each definition can apply to a part of the cell, or the entire cell, beware!
	//This means Properties depend on morphology of the parent cell.
	//optional segment group it applies to, default is the entire cell
	auto aSeg = eAppliedOn.attribute("segment");
	auto aGroup = eAppliedOn.attribute("segmentGroup");
	
	if(aSeg && aGroup){
		log.error(eAppliedOn, "both segment and segmentGroup specified");
		return false;
	}
	else if(aSeg){
		
		Int nml_segid = -1;
		if( !StrToL(aSeg.value(), nml_segid) ){
			log.error(eAppliedOn, "invalid segment id %s", aSeg.value());
			return false;
		}
		//check if it exists in morphology
		Int internal_segid = morph.lookupSegId(nml_segid);
		if( internal_segid < 0 ){
			log.error(eAppliedOn, "segment %ld does not exist in associated Morphology", nml_segid);
			return false;
		}
		
		applied_on.Segment(internal_segid);
		return true;
	}
	else{
		return ParseAcrossSegGroup(log, aGroup.value(), eAppliedOn, morph, applied_on);
	}
}

// get the scaling factor that applies, compared to native units, for that unit name
bool ValidateGetUnits(const ILogProxy &log, const DimensionSet &dimensions, const Dimension &dimension, const char *unit_name, LemsUnit &units){
	if(!dimensions.Has(dimension)){
		log.error("there are no specified %s units", dimensions.Stringify(dimension).c_str() );
		return false;
	}
	for( auto scale : dimensions.GetUnits(dimension) ){
		if(strcmp(unit_name, scale.name.c_str()) == 0){
			units = scale;
			return true;
		}
	}
	//unit name not found in list!
	std::string known_list; for( auto scale : dimensions.GetUnits(dimension) ) known_list += " ", known_list += scale.name;
	log.error("unknown units: %s for %s (supported:%s)", unit_name, dimensions.Stringify(dimension).c_str(), known_list.c_str() );
	return false;
}

//NeuroML physical quantities consist of a numeric, along with an unit name (such as meter, kilometer, etc.) qualifying the quantity the numeric represents. So NeuroML reader code has to check the unit name, to properly read the quantity.

// https://stackoverflow.com/questions/3768862/c-single-template-specialisation-with-multiple-template-parameters
template<typename UnitType, typename ValueType>
struct HelpParseQuantity{ static
bool ParseQuantity(const ImportLogger &log, const pugi::xml_node &eLocation, const char *attr_name, ValueType &num){ 	
	const char *qty_text = eLocation.attribute(attr_name).value();
	double pure_number;
	char unit_name[100]; //if you want more, contact the author for an upgrade
	
	if(!*qty_text){
		log.error(eLocation, "required %s attribute %s missing", UnitType::NAME, attr_name);
		return false;
	}
	
	//number, followed by unit name (which, conveniently, does not comtain spaces)
	if(sscanf(qty_text,"%lf%99s",&pure_number, unit_name) != 2){
		log.error(eLocation, "%s attribute not containing a number and unit", attr_name);
		return false;
	}
	
	//then get the scaling factor that applies, compared to native units, for that unit name
	auto native = Scales<UnitType>::native;
	for( auto scale : Scales<UnitType>::scales ){
		if(strcmp(unit_name, scale.name) == 0){
			
			num = scale.ConvertTo( pure_number, native );
			//printf("valll %f\n",num);
			return true;
		}
	}
	//unit name not found in list!
	std::string known_list;
	for( auto scale : Scales<UnitType>::scales ) known_list += " ", known_list += scale.name;
	log.error(eLocation, "unknown %s attribute units: %s for %s (supported:%s)", attr_name, unit_name, UnitType::NAME, known_list.c_str() );
	return false;
}};
// specialize for unitless, lacking unit markup
template<typename ValueType>
struct HelpParseQuantity<Dimensionless,ValueType>{ static
bool ParseQuantity(const ImportLogger &log, const pugi::xml_node &eLocation, const char *attr_name, ValueType &num){
	
	const char *qty_text = eLocation.attribute(attr_name).value();
	double pure_number;
	
	if(!*qty_text){
		log.error(eLocation, "required %s attribute %s missing", Dimensionless::NAME, attr_name);
		return false;
	}
	
	//just a 𝔭𝔲𝔯𝔢 number
	if(!StrToF(qty_text,pure_number)){
		log.error(eLocation, "%s attribute not containing a \xF0\x9D\x94\xAD\xF0\x9D\x94\xB2\xF0\x9D\x94\xAF\xF0\x9D\x94\xA2 number", attr_name);
		return false;
	}
	
	num = pure_number;
	return true;
}};

template<typename UnitType, typename ValueType>
bool ParseQuantity(const ImportLogger &log, const pugi::xml_node &eLocation, const char *attr_name, ValueType &num){
return HelpParseQuantity<UnitType,ValueType>::ParseQuantity(log, eLocation, attr_name, num); }

template<typename ValueType> 
bool ParseLemsQuantity(const ImportLogger &log, const pugi::xml_node &eLocation, const char *attr_name, const DimensionSet &dimensions, const Dimension dimension, ValueType &num){
	
	const char *qty_text = eLocation.attribute(attr_name).value();
	double pure_number;
	char unit_name[100]; //if you want more, contact the author for an upgrade
	
	if(!*qty_text){
		log.error(eLocation, "required attribute %s missing", attr_name);
		return false;
	}
	
	if(dimension == Dimension::Unity()){
		
		if(!StrToF(qty_text,pure_number)){
			log.error(eLocation, "%s attribute not containing a \xF0\x9D\x94\xAD\xF0\x9D\x94\xB2\xF0\x9D\x94\xAF\xF0\x9D\x94\xA2 number", attr_name);
			return false;
		}
		
		num = pure_number;
		return true;
	}
	
	//number, followed by unit name (which, conveniently, does not comtain spaces)
	if(sscanf(qty_text,"%lf%99s",&pure_number, unit_name) != 2){
		log.error(eLocation, "%s attribute not containing a number and unit", attr_name);
		return false;
	}
	//then get the scaling factor that applies, compared to native units, for that unit name
	LemsUnit scale = dimensions.GetNative(Dimension::Unity());
	if(!ValidateGetUnits(LogWithElement(log,eLocation), dimensions, dimension, unit_name, scale)) return false;
	num = scale.ConvertTo( pure_number, dimensions.GetNative(dimension) );
	//printf("valll %f\n",num);
	return true;
}

bool ParseDimensionAttribute(const ImportLogger &log, const pugi::xml_node &eTag, const DimensionSet &dimensions, Dimension &dimension, const char *sAttributeName = "dimension"){
	auto sDimension = RequiredAttribute(log, eTag, sAttributeName);
	if(!sDimension) return false;
	
	// find the dimension
	if( !dimensions.Has(sDimension) ){
		log.error(eTag, "unknown dimension %s", sDimension);
		return false;
	}
	dimension = dimensions.Get(sDimension);
	return true;
}

template<typename UnitType>
bool ParseValueAcrossSegOrSegGroup(const ImportLogger &log, const pugi::xml_node &eAppliedOn, const char *attr_name, const Morphology &morph , ValueAcrossSegOrSegGroup &applied_on){
	
	if( !ParseQuantity<UnitType>(log, eAppliedOn, attr_name, applied_on.value) ) return false;
	if( !ParseAcrossSegOrSegGroup(log, eAppliedOn, morph, applied_on) ) return false;
	
	return true;
}

bool ParseSpeciesAcrossSegOrSegGroup(const ImportLogger &log, const pugi::xml_node &eSpeciesSpec, const Morphology &morph, const CollectionWithNames<ConcentrationModel> &conc_models,
	CollectionWithNames<IonSpecies> &ion_species, SpeciesAcrossSegOrSegGroup &species_spec){
	
	if( !ParseAcrossSegOrSegGroup( log, eSpeciesSpec, morph, species_spec ) ) return false;
	
	const char *ion_name = eSpeciesSpec.attribute("id").value();
	if( !*ion_name ){
		log.error(eSpeciesSpec, "ion species instance missing id");
		return false;
	}
	const char *ion_name_name = eSpeciesSpec.attribute("ion").value();
	if( *ion_name_name != '\0' && strcmp(ion_name, ion_name_name ) != 0){
		log.error(eSpeciesSpec, "ion species instance id \"%s\" different from ion name \"%s\"",ion_name, ion_name_name);
		return false;
	}
	
	//convert ion species name to id
	species_spec.species = ion_species.idOrNew(ion_name);
	
	if( !ParseQuantity<Concentration>(log, eSpeciesSpec, "initialConcentration", species_spec.initialConcentration) ) return false;
	if( !ParseQuantity<Concentration>(log, eSpeciesSpec, "initialExtConcentration", species_spec.initialExtConcentration) ) return false;
	
	const char *concentration_model_name = eSpeciesSpec.attribute("concentrationModel").value();
	//find the concentration model
	auto conc_model_id = conc_models.get_id(concentration_model_name);
	if(conc_model_id < 0){
		log.error(eSpeciesSpec, "unknown ion concentration model \"%s\"",concentration_model_name);
		return false;
	}
	species_spec.concentrationModel = conc_model_id;
	
	// check if ion types of pool and match too, just for the fun of it 
	if( species_spec.species != conc_models.get(conc_model_id).ion_species ){
		log.error(eSpeciesSpec, "ion species instance different from its concentration model's ion species", concentration_model_name);
		return false;
	}
	
	// LATER also check if a concentration model is already defined, so astonishing things won't happen !
	// silently override concnetration models with last one, for now?
	
	// finished with species specification across part of a cell
	return true;
}

bool ParseBiophysicalProperties(
	const ImportLogger &log, const pugi::xml_node &eBioph, const Morphology &morph, const CollectionWithNames<ConcentrationModel> &conc_models, const CollectionWithNames<IonChannel> &ion_channels,
	const CollectionWithNames<ComponentType> &component_types,
	bool two_ca_pools,
	CollectionWithNames<IonSpecies> &ion_species, BiophysicalProperties &bioph
){
	//parse Morphology and connect it to biophysics
	
	//fill in id, too
	bioph.name = RequiredNmlId(log, eBioph);
	if(!bioph.name) return false;
	
	bioph.two_ca_pools = two_ca_pools;
	
	//required
	const char *membraName = "membraneProperties";
	if(two_ca_pools) membraName = "membraneProperties2CaPools";
	auto eMembr = eBioph.child(membraName);
	if(!eMembr){
		log.error(eBioph, "%s lacks %s", eBioph.name(), membraName);
		return false;
	}
	
	struct ChannelDistributionSubType{
		ChannelDistribution::Type type;
		bool non_uniform;
		bool uses_conductivity;
	};
	//the // marked types are supposed to have varparm elements, in practice only condDensity has ever been used
	const static NameMap<ChannelDistributionSubType> distribution_types = {
		{"channelPopulation"				, {ChannelDistribution::POPULATION	,false, false} }, //
		{"channelDensity"					, {ChannelDistribution::FIXED 		,false, true } }, //
		{"channelDensityVShift"				, {ChannelDistribution::VSHIFT 		,false, true } }, 
		{"channelDensityNernst"				, {ChannelDistribution::NERNST 		,false, true } }, //
		{"ChannelDensityNernstCa2"			, {ChannelDistribution::NERNST_CA2 	,false, true } }, // 
		{"channelDensityGHK"				, {ChannelDistribution::GHK 		,false, false} },
		{"channelDensityGHK2"				, {ChannelDistribution::GHK2 		,false, true } },
		{"channelDensityNonUniform"			, {ChannelDistribution::FIXED 		,true , true } }, //
		{"channelDensityNonUniformNernst"	, {ChannelDistribution::NERNST 		,true , true } }, //
		{"channelDensityNonUniformGHK"		, {ChannelDistribution::GHK 		,true , false} }, //
		// perhaps more in use in the wild ?
	};
	
	
	for(auto eMembrEl: eMembr.children()){
		// first check if it's one of the channel distribution types
		auto disttype_it = distribution_types.find(eMembrEl.name());
		//printf("%s\n",eMembrEl.name());
		if(disttype_it != distribution_types.end()){
			const auto &eDistr = eMembrEl;
			
			ChannelDistribution distribution;
			
			distribution.name = RequiredNmlId(log, eDistr);
			if(!distribution.name) return false;
			
			distribution.type = disttype_it->second.type;
			bool non_uniform = disttype_it->second.non_uniform;
			
			bool uses_conductivity = !(
				distribution.type == ChannelDistribution::POPULATION
				|| distribution.type == ChannelDistribution::GHK
			);
			
			if(uses_conductivity && !non_uniform){
				distribution.conductivity.type = ChannelDistribution::Conductivity::FIXED;
				if( !ParseQuantity<Conductivity>(log, eDistr, "condDensity", distribution.conductivity.value) ) return false;
			}
			
			if(
				distribution.type == ChannelDistribution::FIXED
				|| distribution.type == ChannelDistribution::VSHIFT
				|| distribution.type == ChannelDistribution::POPULATION
			){
				if( !ParseQuantity<Voltage>(log, eDistr, "erev", distribution.erev) ) return false;
			}
			
			if(distribution.type == ChannelDistribution::VSHIFT){
				if( !ParseQuantity<Voltage>(log, eDistr, "vShift", distribution.vshift) ) return false;
			}
			
			if(distribution.type == ChannelDistribution::POPULATION){
				Int population;
				if(!(StrToL(eDistr.attribute("number").value(), population) && population > 0)){
					log.error(eDistr, "ion channel population needs a positive integer \"number\" attribute");
					return false;
				}
				distribution.number = population;
			}
			
			if(distribution.type == ChannelDistribution::GHK){
				if( !ParseQuantity<Permeability>(log, eDistr, "permeability", distribution.permeability) ) return false;
			}
			
			//Common properties are ion and channel
			const char *ion_name = eDistr.attribute("ion").value();
			if( !*ion_name ){
				log.error(eDistr, "ion channel distribution model missing ion species attribute");
				return false;
			}
			
			// if(distribution.type == ChannelDistribution::NERNST && strcmp( eDistr.name(), "ChannelDensityNernstCa2") == 0 ){
			// 	//special case because the user might have input "ca" instead of ca2 ? probably not right
			// 	ion_name = "ca2";
			// }
			
			//convert ion species name to id
			distribution.ion_species = ion_species.idOrNew(ion_name);
			
			
			const char *channel_name = eDistr.attribute("ionChannel").value();
			distribution.ion_channel = ion_channels.get_id(channel_name);
			if(distribution.ion_channel < 0){
				log.error(eDistr, "unknown ion channel type %s", channel_name);
				return false;
			}
			
			//and check with channel's ion too
			if( ion_channels.get(distribution.ion_channel).species != distribution.ion_species
				&& !( (ion_channels.get(distribution.ion_channel).species < 0) && strcmp(ion_name, "non_specific") == 0 ) // in case of non specific
			){
				log.warning(eDistr, "ion type of ion channel distribution %s does not match ion type of ion channel %s", ion_name, channel_name);
				
				// it happens all the time, like when calcium channnels are renamed in distributions,
				// to not affect calcuim-controlled mechanisms
				// proceed 
			}
			
			// and why not check if conductance is defined in case of population distribution
			if( distribution.type == ChannelDistribution::POPULATION ){
				const auto &channel = ion_channels.get(distribution.ion_channel); 
				if( channel.type ==  IonChannel::COMPONENT ){
					const auto &component = component_types.get(channel.component.id_seq);
					if(!( component.exposures.has("g") && component.getExposureDimension("g") == LEMS_Conductance )){
						log.error(eDistr, "ion channel %s does not expose \"g\" of dimension Conductance", channel_name);
						return false;
					}
					log.error(eDistr, "conductance not yet supported for LEMS ion channel %s", channel_name);
					return false;
				}
				else{
					if( std::isnan(channel.conductance) ){
						log.error(eDistr, "ion channel %s has no defined conductance for Population distribution", channel_name);
						return false;
					}
				}
			}
			
			//Now check where the specification applies, it might even be variableParameter!
			if(!non_uniform){
				//the simple case
				
				if(!ParseAcrossSegOrSegGroup(log, eDistr, morph, distribution)) return false;
				
				bioph.membraneProperties.channel_specs.push_back(distribution); //yay!
			}
			else{
				//the split personality case
				// TODO actually why would it have to be a split personality?? check implementation in export, basically if there were multiple distributions with the same name then the lems path would be ambiguous, plus there is no reason that variableparameters would result in multiple distributions.
				for(auto eDistrEl : eDistr.children()){
					if(strcmp(eDistrEl.name(), "variableParameter") == 0){
						const auto &eVarParm = eDistrEl;
						
						// <xs:attribute name="segmentGroup" type="xs:string" use="required"/>
						auto sGroupName = eVarParm.attribute("segmentGroup").value();
						if( !*sGroupName ){
							log.error(eVarParm, "segmentGroup attribute not specified");
							return false;
						}
						if( !ParseAcrossSegGroup(log, sGroupName, eVarParm, morph, distribution) ) return false;
						
						assert( ((const AcrossSegOrSegGroup &) distribution).type == AcrossSegOrSegGroup::Type::GROUP );
						const auto &seg_group = morph.segment_groups[distribution.seqid];
						
						const auto &eInhoVal = eVarParm.child("inhomogeneousValue");
						
						ChannelDistribution::InhomogeneousValue inho;
						
						auto inhoparm_name = eInhoVal.attribute("inhomogeneousParameter").value();
						inho.parm = seg_group.inhomogeneous_parameters.get_id(inhoparm_name);
						if(inho.parm < 0){
							log.error(eInhoVal, "inhomogeneous parameter %s not found in associated segment group", inhoparm_name);
							return false;
						}
						const auto &inhoparm = seg_group.inhomogeneous_parameters.get(inho.parm);
						
						// now validate the formula. Perhaps this will get differentiated according to parameter type LATER
						auto inhoparm_formula = RequiredAttribute(log, eInhoVal, "value");
						if(!inhoparm_formula) return false;
						if( !ParseLemsExpression( inhoparm_formula, inho.value ) ){
							log.error(eInhoVal, "could not parse %s expression", "value");
							return false;
						}
						// Check that the only symbol is the varparm in use
						for( const auto &symname : inho.value.symbol_refs ){
							if( symname != inhoparm.variable ){
								
								log.error(eInhoVal, "unknown expression term %s (the only known one is the variable parameter %s)", 
									symname.c_str(), inhoparm.variable.c_str());
								return false;
							}
						}
						assert( inho.value.symbol_refs.size() <= 1 );
						
						const char *parm = RequiredAttribute(log, eVarParm, "parameter");
						if(!parm) return false;
						
						if( strcmp(parm, "condDensity") == 0 ){
							// see for scaling factor: straight from NEURON 
							// https://github.com/NeuroML/org.neuroml.export/blob/master/src/main/java/org/neuroml/export/neuron/NeuronWriter.java#L1664
							// https://github.com/NeuroML/org.neuroml.export/blob/master/src/main/java/org/neuroml/export/neuron/NeuronWriter.java#L1629
							
							if(!uses_conductivity){
								log.error(eVarParm, " %s cannot accept inhomogeneous %s parameter", eDistr.name(), parm);
								return false;
							}
							
							distribution.conductivity.type = ChannelDistribution::Conductivity::NON_UNIFORM;
							distribution.conductivity.inho = inho;
							// perhaps the value foumula should be parsed here too for consistency LATER?
							// TODO check whether these sub-distributions overlap, and whine!
							// clone this distribution, replacing just the variable fields
							bioph.membraneProperties.channel_specs.push_back(distribution); //ya
						}
						//else if( strcmp(parm, "permeability") == 0 ){
						//	nobody uses anything else but condDensity for now
						//}
						else{
							log.error(eVarParm, "unknown variable parameter %s", parm);
							return false;
						}
						
					}
					else{
						log.error(eDistrEl, "unknown channel density subelement %s", eDistrEl.name());
						return false;
					}
					//y!
				}
			}
		}
		else if(strcmp(eMembrEl.name(), "spikeThresh") == 0){
			ValueAcrossSegOrSegGroup threshold;
			if( !ParseValueAcrossSegOrSegGroup<Voltage>( log, eMembrEl, "value", morph, threshold ) ) return false;
			bioph.membraneProperties.threshold_specs.push_back(threshold);
		}
		else if(strcmp(eMembrEl.name(), "specificCapacitance") == 0){
			ValueAcrossSegOrSegGroup capacitance;
			if( !ParseValueAcrossSegOrSegGroup<SpecificCapacitance>( log, eMembrEl, "value", morph, capacitance ) ) return false;
			bioph.membraneProperties.capacitance_specs.push_back(capacitance);
		}
		else if(strcmp(eMembrEl.name(), "initMembPotential") == 0){
			ValueAcrossSegOrSegGroup initvolt;
			if( !ParseValueAcrossSegOrSegGroup<Voltage>( log, eMembrEl, "value", morph, initvolt ) ) return false;
			bioph.membraneProperties.initvolt_specs.push_back(initvolt);
		}
		else{
			//unknown, ignore
		}
		
	}
	
	//optional
	const char *intracellName = "intracellularProperties";
	if(two_ca_pools) intracellName = "intracellularProperties2CaPools";
	auto eIntra = eBioph.child(intracellName);
	if(eIntra){
		//internal resistivity (ohms per area?) and ionic species present/relevant
		
		//multiple definitions to overlay
		for(auto eIntraEl: eIntra.children()){
			if(strcmp(eIntraEl.name(), "resistivity") == 0){
				ValueAcrossSegOrSegGroup resistivity;
				if( !ParseValueAcrossSegOrSegGroup<Resistivity>( log, eIntraEl, "value", morph, resistivity ) ) return false;
				
				//add to resistivity specifiers for this cell with the known morphology
				bioph.intracellularProperties.resistivity_specs.push_back(resistivity);
			}
			else if(strcmp(eIntraEl.name(), "species") == 0){
				SpeciesAcrossSegOrSegGroup species_spec;
				if(!ParseSpeciesAcrossSegOrSegGroup(log, eIntraEl, morph, conc_models, ion_species, species_spec)) return false;

				bioph.intracellularProperties.ion_species_specs.push_back(species_spec);
			}
			else{
				//unknown, ignore
			}
		}
		
	}
	
	//optional
	auto eExtra = eBioph.child("extracellularProperties");
	if(eExtra){
		for(auto eExtraEl: eExtra.children()){
			if(strcmp(eExtraEl.name(), "species") == 0){
				SpeciesAcrossSegOrSegGroup species_spec;
				if(!ParseSpeciesAcrossSegOrSegGroup(log, eExtraEl, morph, conc_models, ion_species, species_spec)) return false;

				bioph.extracellularProperties.ion_species_specs.push_back(species_spec);
			}
			else{
				//unknown, ignore
				//TODO what about e.g. temperature ?
			}
		}
		
	}
	
	// keep Ca ion, for future reference
	Int Ca_species_seq = ion_species.get_id("ca");
	if( Ca_species_seq < 0 ) Ca_species_seq = ion_species.get_id("Ca");
	bioph.Ca_species_seq = Ca_species_seq;
	
	Int Ca2_species_seq = ion_species.get_id("ca2");
	if( Ca2_species_seq < 0 ) Ca2_species_seq = ion_species.get_id("Ca2");
	bioph.Ca2_species_seq = Ca2_species_seq;
	
	// ignore Standalone child elements
	
	// import complete
	// perform some validation against Morphology:
	//assert(IdListRle::SelfTest());
	IdListRle all_segments = morph.getFullList();
	
	auto Report_MorphoSpecs = [&]( const char * propname, auto specs ){
		IdListRle total_cover;
		
		for(AcrossSegOrSegGroup spec : specs){
			// TODO encapsulate in AcrossSegOrSegGroup utility functions
			if(spec.type == AcrossSegOrSegGroup::SEGMENT) total_cover.Addd(spec.seqid);
			else if(spec.type == AcrossSegOrSegGroup::GROUP) total_cover.Add(morph.segment_groups[spec.seqid].list);
		}
		
		if( total_cover.Count() != (Int) all_segments.Count() ){
			Morphology::SegmentGroup missing;
			missing.list = all_segments.Minus(total_cover);
			
			log.error(eBioph, "%s info is missing for cell segments: %s", propname, morph.Stringify(missing).c_str());
			return false;
		}
		return true;
	};
	
	// Capacitance and initial voltage must have been specified for all segments
	if( !Report_MorphoSpecs("Specific Capacitance", bioph.membraneProperties.capacitance_specs)) return false;
	if( !Report_MorphoSpecs("Initial Voltage", bioph.membraneProperties.initvolt_specs)) return false;
	
	//Resistivity must have been specified all over the cell, for multi-compartment cells
	if(all_segments.Count() > 1){
		if( !Report_MorphoSpecs("Internal Resistivity", bioph.intracellularProperties.resistivity_specs)) return false;
	}
	// XXX check if pools for the same species coincide! What happens then?
	// TODO verify calcium dependencies of ion channels are met, XXX assume zero concentration of missing for now
	
	if(log.debug_printing) bioph.debug_print(morph, ion_species);
	return true;
}

bool ParseQ10( const ImportLogger &log, const pugi::xml_node &eQ10, Q10Settings &q10 ){
	
	
	auto type = RequiredAttribute(log, eQ10, "type");
	if(!type) return false;
	
	if( strcmp(type, "q10Fixed") == 0 ){
		q10.type = Q10Settings::FIXED;
		if(!ParseQuantity<Dimensionless>(log, eQ10, "fixedQ10", q10.q10)) return false;
	}
	else if( strcmp(type, "q10ExpTemp") == 0 ){
		q10.type = Q10Settings::FACTOR;
		if(!ParseQuantity<Dimensionless>(log, eQ10, "q10Factor", q10.q10)) return false;
		if(!ParseQuantity<Temperature>(log, eQ10, "experimentalTemp", q10.experimentalTemp)) return false;
	}
	else{
		log.error(eQ10, "unknown q10Settings type %s", type);
		return false;
	}
	
	return true;
}
	
bool ParseIonChannelBase(const ImportLogger &log, const pugi::xml_node &eChannel, CollectionWithNames<IonSpecies> &ion_species, IonChannel &channel){
	
	auto ion_species_name = eChannel.attribute("species").value();
	if(!*ion_species_name){
		channel.species = -1; // non-specific current
	}
	else{
		channel.species = ion_species.idOrNew(ion_species_name);
	}
	channel.conductance = NAN;
	if(eChannel.attribute("conductance")){
		if(!ParseQuantity<Conductance>(log, eChannel, "conductance", channel.conductance)) return false;
	}
	
	return true;
}

Int ParseComponentInstanceType( const ImportLogger &log, const pugi::xml_node &eInstance, const CollectionWithNames<ComponentType> &component_types, const char *sType = NULL ){
	Int comp_seq = -1;
	if( !*sType || strcmp(sType, "Component" ) == 0 ){
		sType = RequiredComponentType(log, eInstance);
		if(!sType) return -1;
	}
	
	comp_seq = component_types.get_id(sType);
	if( comp_seq < 0 ){
		log.error(eInstance, "unknown component type %s", sType);
		return -1;
	}
	return comp_seq; // yay!
}

bool ParseComponentInstance(const ImportLogger &log, const pugi::xml_node &eInstance, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const char *type, ComponentInstance &instance ){
	
	// just fill in an instance from the given tag, without assuming the component instance is yet ready to instantiate (all constants defined)
	
	auto &comp_seq = instance.id_seq;
	comp_seq = ParseComponentInstanceType(log, eInstance, component_types, type);
	if( comp_seq < 0 ) return false;
	
	const auto &comp = component_types.get(comp_seq);
	
	// check Parameters, for both completeness and possible overrides
	for( auto keyval : comp.properties.names ){
		auto propname = keyval.first;
		const auto &prop = comp.properties.get(keyval.second);
		
		auto sAttr = eInstance.attribute(propname).value();
		if(*sAttr){
			Real override_value;
			if( !ParseLemsQuantity(log, eInstance, propname, dimensions, prop.dimension, override_value ) ) return false;
			instance.parms.push_back( { keyval.second, override_value } );
		}
		
	}
	
	return true;
}

bool ValidateComponentTypeInterface(
	const ImportLogger &log, const pugi::xml_node &eEntityUsing, 
	const ComponentType &comp, const DimensionSet &dimensions, const char *type, 
	const std::map< std::string, ComponentType::Requirement > &provided_requirements, 
	const std::map< std::string, ComponentType::Requirement > &required_exposures,
	const std::map< std::string, ComponentType::EventPortIn > &provided_event_inputs,
	const std::map< std::string, ComponentType::EventPortOut > &required_event_outputs
){
	
	//  each of the component's requirements must be matched by a provided requirement
	for( auto keyval : comp.requirements.names ){
		auto reqname = keyval.first;
		const ComponentType::Requirement &req = comp.requirements.get(keyval.second);
		
		auto it = provided_requirements.find(reqname);
		if( it == provided_requirements.end() ){
			std::string moreinfo;
			moreinfo += "Provided requirements: ";
			for(const auto &keyval : provided_requirements){
				const auto &name = keyval.first;
				const auto &req = keyval.second;
				moreinfo += name + " (" + dimensions.Stringify(req.dimension) + "); ";
			}
			log.error(eEntityUsing, "component %s has requirement %s not satisfied by this component type: %s", type, reqname, moreinfo.c_str());
			return false;
		}
		const ComponentType::Requirement &creq = it->second;
		
		// check dimensions
		if( req.dimension != creq.dimension ){
			log.error(eEntityUsing, "requirement %s dimensions mismatch: expected %s, got %s", reqname, dimensions.Stringify(req.dimension).c_str(), dimensions.Stringify(creq.dimension).c_str());
			return false;
		}
		
		// pass, yay!
	}
	
	for( auto keyval : comp.event_inputs.names ){
		auto inpname = keyval.first;
		// const ComponentType::EventPortIn &evi = comp.event_inputs.get(keyval.second);
		
		auto it = provided_event_inputs.find(inpname);
		if( it == provided_event_inputs.end() ){
			std::string moreinfo;
			moreinfo += "Provided event inputs: ";
			for(const auto &keyval : provided_event_inputs){
				const auto &name = keyval.first;
				// const auto &req = keyval.second;
				moreinfo += name + "; ";
			}
			log.error(eEntityUsing, "component %s has event input %s not satisfied by this component: %s", type, inpname, moreinfo.c_str());
			return false;
		}
		// const ComponentType::EventPortIn &creq = it->second;
		
		// nothing to compare inside them, yay!
		
		// pass, yay!
	}
	
	//  each of the required exposures must be matched by a component's exposure
	for( auto keyval : required_exposures ){
		auto expname = keyval.first.c_str();
		const auto &exp = keyval.second;
		
		auto exp_seq = comp.exposures.get_id(expname);
		if( exp_seq < 0 ){
			log.error(eEntityUsing, "component %s lacks required exposure %s", type, expname);
			return false;
		}
		// const ComponentType::Exposure &cexp = comp.exposures.get(exp_seq);
		
		// check dimensions
		auto edim = comp.getExposureDimension(expname);
		if( exp.dimension != edim ){
			log.error(eEntityUsing, " exposure %s dimensions mismatch: expected %s, got %s", expname, dimensions.Stringify(exp.dimension).c_str(), edim.Stringify().c_str());
			return false;
		}
		// pass, yay!
	}
	
	return true;
}

bool ValidateComponentInstanceCompleteness(const ImportLogger &log, const pugi::xml_node &eInstance, const ComponentType &comp, const char *type, const ComponentInstance &instance){
	
	// check if all parameters are defined, regardless of whether they are ever actually used (they could be used when a conditional fires, but this is not a computable condition to investigate)
	
	std::vector<Real> overridden_parms(comp.properties.contents.size(), NAN);
	for(auto override : instance.parms){
		overridden_parms[override.seq] = override.value;
	}
	
	for( Int seq = 0; seq < (Int)comp.properties.contents.size() ; seq++ ){
		const char *propname = comp.properties.getName(seq);
		const auto &prop = comp.properties.get(seq);
		
		if(!( std::isfinite(prop.value) || std::isfinite(overridden_parms[seq]) ) ){
			log.error(eInstance, "parameter %s of component %s not specified", propname, type);
			return false;
		}
	}
	
	return true;
}

bool ParseInlineComponentInstance(
	const ImportLogger &log, const pugi::xml_node &eInstance, 
	const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, 
	const char *type, 
	const std::map< std::string, ComponentType::Requirement  > &provided_requirements, 
	const std::map< std::string, ComponentType::Requirement  > &required_exposures, 
	const std::map< std::string, ComponentType::EventPortIn  > &provided_event_inputs,
	const std::map< std::string, ComponentType::EventPortOut > &required_event_outputs,
	ComponentInstance &instance 
){
	if( !ParseComponentInstance(log, eInstance, component_types, dimensions, type, instance) ) return false;
	const auto &comp = component_types.get(instance.id_seq);
	if( !ValidateComponentTypeInterface(log, eInstance, comp, dimensions, type, provided_requirements, required_exposures, provided_event_inputs, required_event_outputs ) ) return false;
	if( !ValidateComponentInstanceCompleteness(log, eInstance, comp, type, instance) ) return false;
	
	return true;
}
bool ParseInlineComponentInstance(
	const ImportLogger &log, const pugi::xml_node &eInstance, 
	const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, 
	const char *type, 
	const std::map< std::string, ComponentType::Requirement > &provided_requirements, 
	const std::map< std::string, ComponentType::Requirement > &required_exposures, 
	ComponentInstance &instance 
){
	return ParseInlineComponentInstance(log, eInstance, component_types, dimensions, type, provided_requirements,required_exposures, std::map< std::string, ComponentType::EventPortIn >(), std::map< std::string, ComponentType::EventPortOut >(), instance );
}


std::map< std::string, ComponentType::EventPortOut > required_event_outputs;
// helpers for HH ?
static void CoverCommonRequirement(const char *req_name, Dimension req_dimension, std::map< std::string, ComponentType::Requirement > &provided_requirements){
	ComponentType::Requirement v; v.dimension = req_dimension  ; provided_requirements.insert( std::make_pair(req_name, v) );
};
static void CoverCommonEventIn(const char *req_name, std::map< std::string, ComponentType::EventPortIn > &provided_eventins){
	ComponentType::EventPortIn v; provided_eventins.insert( std::make_pair(req_name, v) );
};
static void CoverCommonEventOut(const char *req_name, std::map< std::string, ComponentType::EventPortOut > &provided_eventouts){
	ComponentType::EventPortOut v; provided_eventouts.insert( std::make_pair(req_name, v) );
};
static void CoverCommonTemperature(std::map< std::string, ComponentType::Requirement > &provided_requirements){
	CoverCommonRequirement( "temperature", LEMS_Temperature, provided_requirements );
}
static void CoverCommonTime(std::map< std::string, ComponentType::Requirement > &provided_requirements){
	CoverCommonRequirement( "t", LEMS_Time, provided_requirements );
}
static void CoverCommonVoltage(std::map< std::string, ComponentType::Requirement > &provided_requirements){
	CoverCommonRequirement( "v", LEMS_Voltage, provided_requirements );
}
static void CoverCommonGlobalRequirements( std::map< std::string, ComponentType::Requirement > &provided_requirements ){
	
	// and don't forget temperature
	CoverCommonTemperature( provided_requirements );
	// and time, why not
	CoverCommonTime( provided_requirements );
	
}
static void CoverCommonIntraCompartmentStuff( std::map< std::string, ComponentType::Requirement > &provided_requirements ){
	
	// first, the globals
	CoverCommonGlobalRequirements( provided_requirements );
	
	CoverCommonVoltage( provided_requirements );
	// since ion channels are specified for multi-segment groups, one could assume calcium concentration exists and resolve it later on;
	// ontherwise, calcium availability could be deduced for that specific group, so error can happen here
	// (but this would cause more confusing error messages in case ion concetration was accidentally omitted for a few segments)
	CoverCommonRequirement( "caConc", LEMS_Concentration, provided_requirements );
	
}
static void CoverCommonRateThingStuff( const IonChannel::Gate::Type gatetype, std::map< std::string, ComponentType::Requirement > &provided_requirements ){
	
	CoverCommonIntraCompartmentStuff( provided_requirements );
	
	if(
		gatetype == IonChannel::Gate::RATES
		|| gatetype == IonChannel::Gate::RATESTAU
		|| gatetype == IonChannel::Gate::RATESTAUINF
	){
		CoverCommonRequirement( "alpha", LEMS_Frequency, provided_requirements );
		CoverCommonRequirement( "beta" , LEMS_Frequency, provided_requirements );
	}
	CoverCommonRequirement ("rateScale", Dimension::Unity(), provided_requirements );
	// TODO if any other case actually happens
}
bool ParseComponentInstanceHHRate(const ImportLogger &log, const pugi::xml_node &eRate, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const IonChannel::Gate::Type gatetype, const char *type, ComponentInstance &instance){
	
	std::map< std::string, ComponentType::Requirement > provided_requirements;
	CoverCommonRateThingStuff( gatetype, provided_requirements );
	
	
	std::map< std::string, ComponentType::Requirement > required_exposures   ;
	ComponentType::Requirement r; r.dimension = LEMS_Frequency; required_exposures   .insert( std::make_pair("r", r) );
	
	if( !ParseInlineComponentInstance(log, eRate, component_types, dimensions, type, provided_requirements, required_exposures, instance) ) return false;
	
	//ComponentType::CommonRequirements gatetime_reqs; // satisfied by core NeuroML
	//gatetime_reqs.membrane_voltage = 1;
	
	/*
	if(!( gatetime_reqs.Superset(comp.common_requirements) )){
		log.error(eRate, "component %s has a requirement not satisfied by this component", type);
		// TODO be more specific
		return false;
	}
	*/
	
	return true;
}

bool ParseComponentInstanceHHTime(const ImportLogger &log, const pugi::xml_node &eRate, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const IonChannel::Gate::Type gatetype, const char *type, ComponentInstance &instance){
	
	std::map< std::string, ComponentType::Requirement > provided_requirements;
	CoverCommonRateThingStuff( gatetype, provided_requirements );
	
	std::map< std::string, ComponentType::Requirement > required_exposures   ;
	ComponentType::Requirement t; t.dimension = LEMS_Time      ; required_exposures   .insert( std::make_pair("t", t) );
	
	if( !ParseInlineComponentInstance(log, eRate, component_types, dimensions, type, provided_requirements, required_exposures, instance) ) return false;
	
	return true;
}

bool ParseComponentInstanceHHVariable(const ImportLogger &log, const pugi::xml_node &eRate, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const IonChannel::Gate::Type gatetype, const char *type, ComponentInstance &instance){
	
	std::map< std::string, ComponentType::Requirement > provided_requirements;
	CoverCommonRateThingStuff( gatetype, provided_requirements );
	
	std::map< std::string, ComponentType::Requirement > required_exposures   ;
	ComponentType::Requirement x; x.dimension = Dimension::Unity() ; required_exposures   .insert( std::make_pair("x", x) );
	
	if( !ParseInlineComponentInstance(log, eRate, component_types, dimensions, type, provided_requirements, required_exposures, instance) ) return false;
	
	return true;
}


bool ParseComponentInstanceConcentrationModel(const ImportLogger &log, const pugi::xml_node &eInstance, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const char *type, ComponentInstance &instance){
	
	std::map< std::string, ComponentType::Requirement > provided_requirements;
	CoverCommonIntraCompartmentStuff( provided_requirements );
	ComponentType::Requirement surfaceArea ; surfaceArea.dimension = LEMS_Area ; provided_requirements.insert( std::make_pair("surfaceArea", surfaceArea) );
	ComponentType::Requirement initialConcentration   ; initialConcentration   .dimension = LEMS_Concentration; provided_requirements.insert( std::make_pair("initialConcentration"   , initialConcentration   ) );
	ComponentType::Requirement initialExtConcentration; initialExtConcentration.dimension = LEMS_Concentration; provided_requirements.insert( std::make_pair("initialExtConcentration", initialExtConcentration) );
	// TODO migrate to an on-demand interface check because ion name may change !
	// or wait for an influx current name common to all ion species! also needed for two ca pools FIXME !!
	ComponentType::Requirement iCa  ; iCa .dimension = LEMS_Current ; provided_requirements.insert( std::make_pair("iCa" , iCa ) );
	ComponentType::Requirement iCa2 ; iCa2.dimension = LEMS_Current ; provided_requirements.insert( std::make_pair("iCa2", iCa2) );
	
	std::map< std::string, ComponentType::Requirement > required_exposures   ;
	ComponentType::Requirement concentration   ; concentration.dimension    = LEMS_Concentration; required_exposures.insert( std::make_pair("concentration"   , concentration   ) );
	ComponentType::Requirement extConcentration; extConcentration.dimension = LEMS_Concentration; required_exposures.insert( std::make_pair("extConcentration", extConcentration) );
	
	
	if( !ParseInlineComponentInstance(log, eInstance, component_types, dimensions, type, provided_requirements, required_exposures, instance) ) return false;
	
	return true;
}

bool ParseComponentInstanceIonChannel(const ImportLogger &log, const pugi::xml_node &eThing, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const char *type, ComponentInstance &instance){
	
	std::map< std::string, ComponentType::Requirement > provided_requirements;
	CoverCommonIntraCompartmentStuff( provided_requirements );
	
	std::map< std::string, ComponentType::Requirement > required_exposures   ;
	// ComponentType::Requirement g; g.dimension = LEMS_Conductance ; required_exposures.insert( std::make_pair("g", g) ); LATER as needed
	ComponentType::Requirement fopen; fopen.dimension = Dimension::Unity() ; required_exposures   .insert( std::make_pair("fopen", fopen) );
	
	if( !ParseInlineComponentInstance(log, eThing, component_types, dimensions, type, provided_requirements, required_exposures, instance) ) return false;
	
	return true;
}
bool ParseComponentInstanceGate(const ImportLogger &log, const pugi::xml_node &eGate, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const char *type, ComponentInstance &instance){
	
	std::map< std::string, ComponentType::Requirement > provided_requirements;
	CoverCommonIntraCompartmentStuff( provided_requirements );
	
	std::map< std::string, ComponentType::Requirement > required_exposures   ;
	ComponentType::Requirement q; q.dimension = Dimension::Unity() ; required_exposures   .insert( std::make_pair("q", q) );
	ComponentType::Requirement fcond; fcond.dimension = Dimension::Unity() ; required_exposures   .insert( std::make_pair("fcond", fcond) );
	
	if( !ParseInlineComponentInstance(log, eGate, component_types, dimensions, type, provided_requirements, required_exposures, instance) ) return false;
	
	return true;
}
bool ParseComponentInstanceConductanceScaling(const ImportLogger &log, const pugi::xml_node &eInstance, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const char *type, ComponentInstance &instance){
	
	std::map< std::string, ComponentType::Requirement > provided_requirements;
	CoverCommonIntraCompartmentStuff( provided_requirements );
	
	std::map< std::string, ComponentType::Requirement > required_exposures   ;
	ComponentType::Requirement factor; factor.dimension = Dimension::Unity() ; required_exposures.insert( std::make_pair("factor", factor) );
	
	if( !ParseInlineComponentInstance(log, eInstance, component_types, dimensions, type, provided_requirements, required_exposures, instance) ) return false;
	
	return true;
}
bool ParseComponentInstanceInputCurrent(const ImportLogger &log, const pugi::xml_node &eInstance, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const char *type, ComponentInstance &instance){
	
	std::map< std::string, ComponentType::Requirement > provided_requirements;
	CoverCommonIntraCompartmentStuff( provided_requirements );
	
	std::map< std::string, ComponentType::Requirement > required_exposures   ;
	// all required exposures are determined by interface matching, at instantiation time
	
	if( !ParseInlineComponentInstance(log, eInstance, component_types, dimensions, type, provided_requirements, required_exposures, instance) ) return false;
	
	return true;
}

bool ParseComponentInstanceArtificialCell(const ImportLogger &log, const pugi::xml_node &eInstance, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const char *type, ComponentInstance &instance){
	
	std::map< std::string, ComponentType::Requirement > provided_requirements;
	CoverCommonGlobalRequirements( provided_requirements );
	// TODO perhaps be pickier and resolve component type, to determine additional specializations?
	// if( !ParseComponentInstance(log, eInstance, component_types, dimensions, type, instance) ) return false;
	CoverCommonRequirement( "iSyn", LEMS_Current, provided_requirements );
	CoverCommonRequirement( "ISyn", Dimension::Unity(), provided_requirements );
	
	std::map< std::string, ComponentType::Requirement > required_exposures   ;
	// situational requirements will be resolved at synaptic/input projection time
	
	std::map< std::string, ComponentType::EventPortIn > provided_event_inputs;
	std::map< std::string, ComponentType::EventPortOut > required_event_outputs;
	// situational requirements will be resolved at synaptic/input projection time
	
	CoverCommonEventIn( "in", provided_event_inputs ); // just in case, no harm in adding a placeholder
	
	
	if( !ParseInlineComponentInstance(log, eInstance, component_types, dimensions, type, provided_requirements, required_exposures, provided_event_inputs, required_event_outputs, instance) ) return false;
	
	return true;
}
bool ParseComponentInstanceSynapticComponent(const ImportLogger &log, const pugi::xml_node &eInstance, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const char *type, ComponentInstance &instance){
	
	std::map< std::string, ComponentType::Requirement > provided_requirements;
	CoverCommonIntraCompartmentStuff( provided_requirements );
	
	std::map< std::string, ComponentType::Requirement > required_exposures   ;
	CoverCommonRequirement( "i", LEMS_Current, required_exposures );
	
	std::map< std::string, ComponentType::EventPortIn > provided_event_inputs;
	std::map< std::string, ComponentType::EventPortOut > required_event_outputs;
	
	// TODO allow synapses exposing only conductance & erev to pass
	
	// NB need to resolve component type, to determine additional specializations
	// (will be redundantly resolved in ParseInlineComponentInstance)
	if( !ParseComponentInstance(log, eInstance, component_types, dimensions, type, instance) ) return false;
	
	// switch whether vpeer or spike in is provided ? i guess 
	// const auto &comptype = component_types.get(instance.id_seq);
	// TODO actually detect if there is an inconsistency here
	// see what to do about interface matching, probably put the hasvpeer/spikein checks there at syn.connection time
	//if( comptype.HasVpeer(component_types) ){
		CoverCommonRequirement( "vpeer", LEMS_Voltage, provided_requirements );
	//}
	
	//if( comptype.HasSpikeIn(component_types) ){
		CoverCommonEventIn( "in", provided_event_inputs );
	//}
	
	if( !ParseInlineComponentInstance(log, eInstance, component_types, dimensions, type, provided_requirements, required_exposures, provided_event_inputs, required_event_outputs, instance) ) return false;
	
	return true;
}

bool parseHHVariable(const ImportLogger &log, const pugi::xml_node &eRate, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const IonChannel::Gate::Type gatetype, IonChannel::Rate &rate){
	
	auto type = eRate.attribute("type").value();
	if(!*type){
		log.error(eRate, "variable formula requires type attribute");
		return false;
	}
	
	const static NameMap<IonChannel::Rate::Type> gate_types = {
		{"HHExpVariable"		, IonChannel::Rate::EXPONENTIAL	},
		{"HHExpLinearVariable"	, IonChannel::Rate::EXPLINEAR	},
		{"HHSigmoidVariable"	, IonChannel::Rate::SIGMOID		},
	};
	auto it = gate_types.find(type);
	
	if(it == gate_types.end()){
		rate.type = IonChannel::Rate::COMPONENT;
		
		return ParseComponentInstanceHHVariable(log, eRate, component_types, dimensions, gatetype, type, rate.component);
	}
	else{
		rate.type = it->second;
		if(!ParseQuantity<Dimensionless>(log, eRate, "rate", rate.formula.rate)) return false;
		if(!ParseQuantity<Voltage>(log, eRate, "midpoint", rate.formula.midpoint)) return false;
		if(!ParseQuantity<Voltage>(log, eRate, "scale", rate.formula.scale)) return false;
		return true;
	}
	
}
bool parseHHRate(const ImportLogger &log, const pugi::xml_node &eRate, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const IonChannel::Gate::Type gatetype, IonChannel::Rate &rate){
	
	auto type = eRate.attribute("type").value();
	if(!*type){
		log.error(eRate, "rate requires type attribute");
		return false;
	}
	
	const static NameMap<IonChannel::Rate::Type> gate_types = {
		{"HHExpRate"		, IonChannel::Rate::EXPONENTIAL	},
		{"HHExpLinearRate"	, IonChannel::Rate::EXPLINEAR		},
		{"HHSigmoidRate"	, IonChannel::Rate::SIGMOID		},
	};
	auto it = gate_types.find(type);
	
	if(it == gate_types.end()){
		rate.type = IonChannel::Rate::COMPONENT;
		
		return ParseComponentInstanceHHRate(log, eRate, component_types, dimensions, gatetype, type, rate.component);
	}
	else{
		rate.type = it->second;
		
		if(!ParseQuantity<Frequency>(log, eRate, "rate", rate.formula.rate)) return false;
		if(!ParseQuantity<Voltage>(log, eRate, "midpoint", rate.formula.midpoint)) return false;
		if(!ParseQuantity<Voltage>(log, eRate, "scale", rate.formula.scale)) return false;
		return true;
	}
	
}
bool parseHHTime(const ImportLogger &log, const pugi::xml_node &eRate, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const IonChannel::Gate::Type gatetype, IonChannel::Rate &rate){
	
	auto type = eRate.attribute("type").value();
	if(!*type){
		log.error(eRate, "rate requires type attribute");
		return false;
	}
	
	const static NameMap<IonChannel::Rate::Type> gate_types = {
		{"fixedTimeCourse"	, IonChannel::Rate::FIXED	},
	};
	auto it = gate_types.find(type);
	
	if(it == gate_types.end()){
		rate.type = IonChannel::Rate::COMPONENT;
		
		return ParseComponentInstanceHHTime(log, eRate, component_types, dimensions, gatetype, type, rate.component);
	}
	else{
		rate.type = it->second;
		if(!ParseQuantity<Time>(log, eRate, "tau", rate.formula.constant)) return false;
		return true;
	}
	
}

//helpers for when they come in pairs
bool parseHHTauInf(const ImportLogger &log, const pugi::xml_node &eWithRates, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const IonChannel::Gate::Type gatetype, IonChannel::Rate &tau, IonChannel::Rate &inf){
	
	auto eInf = eWithRates.child("steadyState"), eTau = eWithRates.child("timeCourse");
	if(!( eInf && eTau )){
		log.error(eWithRates, "must have steadyState and timeCourse");
		return false;
	}
	if(!parseHHVariable(log, eInf, component_types, dimensions, gatetype, inf)) return false;
	if(!parseHHTime(log, eTau, component_types, dimensions, gatetype, tau)) return false;
	
	return true;
}

bool parseHHForRev(const ImportLogger &log, const pugi::xml_node &eWithRates, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, const IonChannel::Gate::Type gatetype, IonChannel::Rate &a, IonChannel::Rate &b){
	
	auto eAlpha = eWithRates.child("forwardRate"), eBeta = eWithRates.child("reverseRate");
	if(!( eAlpha && eBeta )){
		log.error(eWithRates, "must have forwardRate and reverseRate");
		return false;
	}
	if(!parseHHRate(log, eAlpha, component_types, dimensions, gatetype, a)) return false;
	if(!parseHHRate(log, eBeta,  component_types, dimensions, gatetype, b)) return false;
	
	return true;
}

bool parseKSTransitionBase(const ImportLogger &log, const pugi::xml_node &eTransition, const IonChannel::GateKS &ks, IonChannel::GateKS::Transition_Base &fromto){
	
	
	auto sFrom = eTransition.attribute("from").value(), sTo = eTransition.attribute("to").value();
	if(!( *sFrom && *sTo )){
		log.error(eTransition, "transition requires from, to attributes");
		return false;
	}
	if(!ks.state_names.count(sFrom)){
		log.error(eTransition, "\'from\' state %s not found", sFrom);
		return false;
	}
	if(!ks.state_names.count(sTo)){
		log.error(eTransition, "\'to\' state %s not found", sTo);
		return false;
	}
	fromto.from = ks.state_names.find(sFrom)->second;
	fromto.to   = ks.state_names.find(sTo  )->second;
	
	return true;
}

bool ParseIonChannel(const ImportLogger &log, const pugi::xml_node &eChannel, const CollectionWithNames<ComponentType> &component_types, const DimensionSet &dimensions, CollectionWithNames<IonSpecies> &ion_species, IonChannel &channel){
	
	if(!ParseIonChannelBase(log, eChannel, ion_species, channel)) return false;
	
	// determine type of ion channel's children
	// (could they ever be different types each? only in-the-wild models may tell)
	// for now, yell if they are not the same type
	
	const static NameMap<IonChannel::Gate::Type> gate_types = {
		{"gate"						, IonChannel::Gate::NONE			},
		{"gateHHrates"				, IonChannel::Gate::RATES			},
		{"gateHHtauInf"				, IonChannel::Gate::TAUINF			},
		{"gateHHratesTau"			, IonChannel::Gate::RATESTAU		},
		{"gateHHratesInf"			, IonChannel::Gate::RATESINF		},
		{"gateHHratesTauInf"		, IonChannel::Gate::RATESTAUINF		},
		{"gateHHInstantaneous"		, IonChannel::Gate::INSTANTANEOUS	},
		{"gateFractional"			, IonChannel::Gate::FRACTIONAL		},
		// just admit you want whatever type of gate, wherever
		// why not put a tree structure of sub-gates while you're at it
		// perhaps you should leave such horrors for LEMS? please?
		// whatever, support it since undetermined gate implies it occurring
		// and there is no substantial reason it could not be here
		{"gateKS"					, IonChannel::Gate::KINETIC			},
	};
	
	// channel.type = IonChannel::PASSIVE; // implicitly
	// bool verified_passive = false; // or explicitly
	// neither, it's just an ion channel with no gates at all
	
	// get some type information	
	const char *chantype_str = TagNameOrType(eChannel);
	
	// check if it is a lems type, merge with inner core NeuroML elements MUCH LATER
	if( 
		component_types.has( chantype_str )
		&& component_types.get( chantype_str ).extends == ComponentType::ION_CHANNEL 
	){
		channel.type = IonChannel::COMPONENT;
		
		if( !ParseComponentInstanceIonChannel(log, eChannel, component_types, dimensions, chantype_str, channel.component) ) return false;
		
		return true; // yay!
	}
	
	// otherwise it is a native type
	channel.type = IonChannel::NATIVE;
	
	/*
	if( *chantype_str){
		if(strcmp(chantype_str, "ionChannelHH")){
			//of course it is
		}
		else if(strcmp(chantype_str, "ionChannelPassive")){
			verified_passive = true; // now we know for sure, clash if gates appears further on
		}
		else{
			log.error(eChannel, "unknown ion channel type %s", chantype_str);
			return false;
		}
	}
	*/
	for(auto eChanElm : eChannel.children()){
		
		const char *sType = TagNameOrType(eChanElm);
		
		bool is_gate = false;
		bool is_scaling = false;
		
		// TODO check this too
		// https://github.com/OpenSourceBrain/NSGPortalShowcase/blob/9e8b73cd18854e5d530b04b884603d076b661b6e/examples/NetPyNE2/TestL23/TestL23/kc.channel.nml#L45
		
		// what is this child element ?
		// perhaps annotation stuff?
		// NOTE wherever there is Standalone, notes, annotation and property tage may exist. Handle properties as they appear, for they must be documented to be used.
		if(
			strcmp(eChanElm.name(), "notes") == 0
			|| strcmp(eChanElm.name(), "annotation") == 0
		){
			continue;
		}
		// is it a scaling factor?
		if( strcmp( sType, "q10ConductanceScaling" ) == 0 ){
			is_scaling = true;
		}
		else if( 
			component_types.has(sType)
			&& component_types.get(sType).extends == ComponentType::CONDUCTANCE_SCALING 
		){
			is_scaling = true;
		}
		
		// is it a gate?
		if( gate_types.count(sType)  > 0 ){
			is_gate = true;
		}
		else if( 
			component_types.has(sType) 
			&& component_types.get(sType).extends == ComponentType::GATE 
		){
			is_gate = true;
		}
		if( is_scaling ){
			const auto &eScaling = eChanElm;
			
			if( channel.conductance_scaling.type != IonChannel::ConductanceScaling::NONE ){
				log.error(eScaling, "conductance scaling already defined once");
				return false;
			}
			
			if( strcmp( sType, "q10ConductanceScaling" ) == 0 ){
				channel.conductance_scaling.type = IonChannel::ConductanceScaling::Q10;
				
				auto &q10 = channel.conductance_scaling.q10;
				q10.type = Q10Settings::FACTOR;
				if(!ParseQuantity<Dimensionless>(log, eScaling, "q10Factor", q10.q10)) return false;
				if(!ParseQuantity<Temperature>(log, eScaling, "experimentalTemp", q10.experimentalTemp)) return false;
				
			}
			else{
				// should be a LEMS component
				channel.conductance_scaling.type = IonChannel::ConductanceScaling::COMPONENT;
				if ( !ParseComponentInstanceConductanceScaling( log, eScaling, component_types, dimensions, sType, channel.conductance_scaling.component ) ) return false;
				
			}
			
			// channel conductance scaling has been set, yay!
			continue;
		}
		if( is_gate ){
			const auto &eGate = eChanElm;
			IonChannel::Gate gate;
			
			auto &gatetype = gate.type;
			gatetype = IonChannel::Gate::NONE;
			
			const char *gatetype_str = eGate.name();
			
			auto it = gate_types.find( gatetype_str );
			if(it == gate_types.end()){
				
				// perhaps LEMS component? 
				gatetype = IonChannel::Gate::NONE;
			}
			else{
				gatetype = it->second;
			}
			
			
			
			auto sGateId = RequiredNmlId(log, eGate);
			if(!sGateId) return false;
			
			
			if( gatetype == IonChannel::Gate::NONE ){
				
				// if it is "gate", must check "type" attribute
				if( strcmp(gatetype_str, "gate") == 0){
					gatetype_str = eGate.attribute("type").value();
					if( !*gatetype_str ){
						log.error(eGate, "must specify \"type\" attribute");
						return false;
					}
				}
				else{
					// it is what the tag name is
				}
				
				auto itt = gate_types.find(gatetype_str);
				if(itt == gate_types.end()){
					if(strcmp(gatetype_str, "gateKS") == 0){
						
						//log.error(eGate, "please use ionChannelKS (or LEMS) for kinetic-scheme gates");
						//return false;
						
						//nah, fall through without complaining
					}
					else{
						// perhaps a LEMS component then?
						gatetype = IonChannel::Gate::COMPONENT;
						if ( !ParseComponentInstanceGate( log, eGate, component_types, dimensions, gatetype_str, gate.component ) ) return false;
						
						channel.gates.add(gate, sGateId);
						continue;
					}
				}
				gatetype = itt->second;
				if(gatetype == IonChannel::Gate::NONE){
					log.error(eGate, "gate type \"%s\" is too abstract", gatetype_str);
					return false;
				}
			}
			
			Int instances;
			if(!RequiredPositiveAttribute(log, eGate, "instances", instances)) return false;
			
			// parse and add gate! yay!
			if(gatetype == IonChannel::Gate::INSTANTANEOUS){
				// just its instantaneous value
				gate.instantaneous.instances = instances;
				auto eSS = eGate.child("steadyState");
				if(!eSS){
					log.error(eGate, "instantaneous gate lacks steadyState");
					return false;
				}
				if(!parseHHVariable(log, eSS, component_types, dimensions, gatetype, gate.instantaneous.steadyState)) return false;
				
				channel.gates.add(gate, sGateId);
				continue;
			}
			
			// check Q10
			auto ParseQ10Optional = []( const ImportLogger &log, const pugi::xml_node &eTempDependent, Q10Settings &q10 ){
	
				auto eQ10 = eTempDependent.child("q10Settings");
				if( eQ10 ){
					return ParseQ10( log, eQ10, q10 );
				}
				else{
					// default
					q10.type = Q10Settings::FIXED;
					q10.q10 = 1;
					
					//log.error(eTempDependent, "q10Settings element is required");
					// return false;
				}
				
				return true;
			};
			
			Q10Settings q10;
			if( !ParseQ10Optional(log, eGate, q10) ) return false;
			
			if(gatetype == IonChannel::Gate::FRACTIONAL){
				IonChannel::GateFractional fga;
				fga.instances = instances;
				fga.q10 = q10;
				
				for( auto eGateEl : eGate.children()){
					if( strcmp( eGateEl.name(), "subGate") == 0 ){
						IonChannel::SubGateFractional sga;
						if( !ParseQ10Optional(log, eGateEl, sga.q10) ) return false;
						if( !ParseQuantity<Dimensionless>(log, eGateEl, "fractionalConductance", sga.fraction_of_conductivity) ) return false;
						if( !parseHHTauInf(log, eGateEl, component_types, dimensions, gatetype, sga.timeCourse, sga.steadyState) ) return false;
						fga.subgates.push_back(sga);
					}
					else{
						// unknown, let it pass
					}
				}
				
				//sanity check
				if(fga.subgates.empty()){
					log.error(eGate, "must have at least one fractional gate");
					return false;
				}
				
				//finished, yay!
				gate.fractional = channel.fractional_gates.size();
				channel.fractional_gates.push_back(fga);
				channel.gates.add(gate, sGateId);
				continue;
			}
			
			if(gatetype == IonChannel::Gate::KINETIC){
				IonChannel::GateKS ks;
				ks.instances = instances;
				ks.q10 = q10;
				
				// since the NeuroML authors wanted to be mean and split forward from reverse transitions
				// and the NEURON exporter does not trust the user-provided ordering, I should ensure proper ordering too
				
				struct ForRevTransitionRef{
					const char *from, *to;
					pugi::xml_node node;
					bool operator<( const ForRevTransitionRef &rhs) const {
						auto comp1 = strcmp(from, rhs.from);
						if(comp1) return (comp1 < 0);
						else return ( strcmp(to, rhs.to) < 0);
					};
					bool operator==( const ForRevTransitionRef &rhs) const {
						return strcmp(from, rhs.from) == 0
							&& strcmp(to, rhs.to) == 0;
					}
				};
				std::vector<ForRevTransitionRef> forward_transitions, reverse_transitions;
				
				for( auto eGateEl : eGate.children()){
					
					if( strcmp( eGateEl.name(), "closedState") == 0 || strcmp( eGateEl.name(), "openState") == 0){
						auto name = RequiredNmlId(log, eGateEl);
						if(!name) return false;
						if(ks.state_names.count(name) > 0){
							log.error(eGateEl, "state %s already defined", name);
							return false;
						}
						
						Int new_id = ks.state_names.size();
						ks.state_names.insert(std::make_pair(name, new_id));
						if(strcmp( eGateEl.name(), "closedState") == 0){
							ks.closed_states.push_back(new_id);
						}
						else{
							ks.open_states.push_back(new_id);
						}
					}
					else if( strcmp( eGateEl.name(), "forwardTransition") == 0 ){
						ForRevTransitionRef ref = { eGateEl.attribute("from").value(), eGateEl.attribute("to").value(), eGateEl };
						forward_transitions.push_back(ref);// handle it later
					}
					else if( strcmp( eGateEl.name(), "reverseTransition") == 0 ){
						ForRevTransitionRef ref = { eGateEl.attribute("from").value(), eGateEl.attribute("to").value(), eGateEl };
						reverse_transitions.push_back(ref);// handle it later
					}
					else if( strcmp( eGateEl.name(), "tauInfTransition") == 0 ){
						IonChannel::GateKS::Transition tt;
						tt.type = IonChannel::GateKS::Transition::TAU_INF;
						
						if( !parseKSTransitionBase(log, eGateEl, ks, tt.tauinf) ) return false;
						if( !parseHHTauInf(log, eGateEl, component_types, dimensions, gatetype, tt.tauinf.steadyState, tt.tauinf.timeCourse) ) return false;
						
						ks.transitions.push_back(tt); //yay!
					}
					else{
						// unknown, let it pass
					}
				}
				
				//sanity check
				if( ks.open_states.empty() && ks.closed_states.empty() ){
					log.error(eGate, "must have at least one state");
					return false;
				}
				
				//now order the transitions so forward and reverse are side-by-side in the arrays
				std::sort(forward_transitions.begin(), forward_transitions.end());
				std::sort(reverse_transitions.begin(), reverse_transitions.end());
				
				//check for duplicates too
				for(size_t i = 1; i < forward_transitions.size(); i++){
					if(forward_transitions[i - 1] == forward_transitions[i]){
						log.error(forward_transitions[i-1].node, "duplicate definition for %s->%s transition", forward_transitions[i].from, forward_transitions[i].to);
						return false;
					}
				}
				for(size_t i = 1; i < reverse_transitions.size(); i++){
					if(reverse_transitions[i - 1] == reverse_transitions[i]){
						log.error(reverse_transitions[i-1].node, "duplicate definition for %s<-%s transition", reverse_transitions[i].from, reverse_transitions[i].to);
						return false;
					}
				}
				
				if(forward_transitions.size() != reverse_transitions.size()){
					log.error(eGate, "forward transitions are %zd while reverse transitions are %zd", forward_transitions.size(), reverse_transitions.size());
					return false;
				}
				
				for(size_t i = 0; i < forward_transitions.size(); i++){
					const ForRevTransitionRef &f = forward_transitions[i], &r = reverse_transitions[i];
					// printf(" %s->%s vs %s<-%s\n", f.from, f.to, r.from, r.to);
					
					if(!( f == r )){
						if(f < r){
							log.error(eGate, "missing reverse transition for %s<-%s", r.from, r.to);
						}
						else{
							log.error(eGate, "missing forward transition for %s->%s", f.from, f.to);
						}
						return false;
					}
					
					IonChannel::GateKS::Transition tt;
					tt.type = IonChannel::GateKS::Transition::FORWARD_REVERSE;
					
					if( !parseKSTransitionBase(log, forward_transitions[i].node, ks, tt.forrev) ) return false;
					auto eAlpha = forward_transitions[i].node.child("rate");
					if(!eAlpha){
						log.error(forward_transitions[i].node, "requires <rate> element");
						return false;
					}
					auto eBeta = reverse_transitions[i].node.child("rate");
					if(!eBeta){
						log.error(reverse_transitions[i].node, "requires <rate> element");
						return false;
					}
					if( !parseHHRate(log, eAlpha, component_types, dimensions, gatetype, tt.forrev.forwardRate) ) return false;
					if( !parseHHRate(log, eBeta , component_types, dimensions, gatetype, tt.forrev.reverseRate) ) return false;
					
					ks.transitions.push_back(tt); //yay!
				}
				
				//finished, yay!
				gate.kinetic = channel.kinetic_gates.size();
				channel.kinetic_gates.push_back(ks);
				channel.gates.add(gate, sGateId);
				continue;
			}
			
			//otherwise it is some alpha-beta-tau-inf mix
			gate.gaga.instances = instances;
			gate.gaga.q10 = q10;
			
			if(
				gatetype == IonChannel::Gate::TAUINF
				|| gatetype == IonChannel::Gate::RATESTAU
				|| gatetype == IonChannel::Gate::RATESTAUINF
			){
				auto eTau = eGate.child("timeCourse");
				if(!( eTau )){
					log.error(eGate, "must have timeCourse");
					return false;
				}
				if(!parseHHTime(log, eTau, component_types, dimensions, gatetype, gate.gaga.timeCourse)) return false;
			}
			
			if(
				gatetype == IonChannel::Gate::TAUINF
				|| gatetype == IonChannel::Gate::RATESINF
				|| gatetype == IonChannel::Gate::RATESTAUINF
			){
				auto eInf = eGate.child("steadyState");
				if(!( eInf )){
					log.error(eGate, "must have steadyState");
					return false;
				}
				if(!parseHHVariable(log, eInf, component_types, dimensions, gatetype, gate.gaga.steadyState)) return false;
			}
			
			if(
				gatetype == IonChannel::Gate::RATES
				|| gatetype == IonChannel::Gate::RATESTAU
				|| gatetype == IonChannel::Gate::RATESINF
				|| gatetype == IonChannel::Gate::RATESTAUINF
			){
				if( !parseHHForRev(log, eGate, component_types, dimensions, gatetype, gate.gaga.forwardRate, gate.gaga.reverseRate) ) return false;
			}
			
			channel.gates.add(gate, sGateId);
			continue;
		}
		
	}
	
	return true;
}

// merge this ref parsing with LemaQuantityPath processing LATER
// Fails without writing to log, just returns false
bool ParseSynapseCellRef(const char *refspec, Int &id){
	
	// Getting the cell reference is tons of fun: https://github.com/NeuroML/NeuroML2/issues/109
	// Use the official parsing algorithm, which just happens to be the jLEMS source
	// https://github.com/NeuroML/org.neuroml.export/blob/development/src/main/java/org/neuroml/export/utils/Utils.java#L119
	
	const char *ptr = refspec; // indexing progress through parsing the string
	//printf("start %s\n",ptr);
	if(strncmp(refspec, "../", 3) == 0){
		
		ptr += 3; // skip past ../
		//printf("%s\n",ptr);
		// ptr = strchr(ptr, '/'); //and what follows
		// if(!ptr) return false;
		// ptr++; //skip past / found
	}
	//printf("absolute %s\n",ptr);
	// now ptr is either of the form pop/<number>/ ... or the form pop[number]/... or in the new format, just the instance <number> itself
	//skip past / or [
	auto bracket = strchr(ptr, '[');
	auto slash =  strchr(ptr,'/');
	if(bracket){
		ptr = bracket + 1;
	}
	else if(slash){
		ptr = slash + 1;
	}
	//printf("position %s\n",ptr);
	return StrToL(ptr, id, false);
}
bool ParseSynapsePopulationRef(const char *refspec, std::string &pop_name){
	const char *ptr = refspec; // indexing progress through parsing the string
	const char *endptr = refspec + strlen(refspec);
	//printf("start %s\n",ptr);
	if(strncmp(refspec, "../", 3) == 0){
		ptr += 3; // skip past ../
	}
	auto bracket = strchr(ptr, '[');
	if(bracket && bracket < endptr) endptr = bracket;
	auto slash = strchr(ptr, '/');
	if(slash   && slash   < endptr) endptr = slash  ;
	
	pop_name = std::string(ptr, endptr - ptr);
	return true;
}

bool parseCompartmentTarget(const ImportLogger &log, const pugi::xml_node &eTarget,
	const CollectionWithNames<Morphology> &morphologies,
	const CollectionWithNames<CellType> &cell_types,
	const CollectionWithNames<Network::Population> &populations,
	Int &population_seq, Int &instance_seq, Int &segment_seq, Real &fractionAlong,
	bool known_population = false // e.g. in cases like inputList
){
	auto targetPath = RequiredAttribute(log, eTarget, "target");
	if(!targetPath) return false;
	
	std::string pop_name;
	Int input_cell_instance_id;
	if(!(
		ParseSynapsePopulationRef(targetPath, pop_name)
		&& ParseSynapseCellRef(targetPath, input_cell_instance_id)
	)){
		log.error(eTarget, "target path %s not valid", targetPath);
		return false;
	}
	
	// and validate
	if(known_population){
		// no need to find population, it might even be missing
	}
	else{
		population_seq = populations.get_id(pop_name.c_str());
		if(population_seq < 0){
			log.error(eTarget, "target population %s not found", pop_name.c_str());
			return false;
		}
	}
	// TODO: validate the name of the population in any case, because if there is a discrepancy
	// jNML will prefer the population name on the XML path and thus give a different result !!
	// DataWritper paths are already checked for this iirc
	
	const Network::Population &incident_population = populations.get(population_seq);
	
	instance_seq = incident_population.instances.getSequential(input_cell_instance_id); //convert to sequential form
	if(instance_seq < 0){
		log.error(eTarget, "cell instance id %ld not present in target population", input_cell_instance_id);
		return false;
	}
	
	// also get segment id
	Int segId = 0;
	auto *sIncident = eTarget.attribute("segmentId").value();
	if( *sIncident ){
		if( !StrToL(sIncident, segId) ){
			log.error(eTarget, "target segment id must be a positive integer");
			return false;
		}	
	}
	
	auto cell_type = cell_types.get( incident_population.component_cell );
	if( cell_type.type == CellType::PHYSICAL ){
		const Morphology &incident_morph = morphologies.get( cell_type.physical.morphology );
		segment_seq = incident_morph.segments.getSequential(segId); //convert to sequential form
	}
	else if( cell_type.type == CellType::ARTIFICIAL ){
		// point neurons, they have only one segment, with id 0
		if( segId == 0 ) segment_seq = 0;
		else return false;// segment_seq = -1;
	}
	else{
		log.error(eTarget, "internal error: compartment target: cell type %d", cell_type.type);
		return false;
	}
	if(segment_seq < 0){
		log.error(eTarget, "%s segment id %ld not present in target cell %ld", (*sIncident)?"implied ":"", segId, input_cell_instance_id);
		return false;
	}
	
	// and fractionAlong
	if(eTarget.attribute("fractionAlong")){
		if( !ParseQuantity<Dimensionless>(log, eTarget, "fractionAlong", fractionAlong) ) return false;
		if(!( 0 <= fractionAlong && fractionAlong <= 1.0 )){
			log.error(eTarget, "fractionAlong not between 0 and 1", segId, input_cell_instance_id);
			return false;
		}
	}
	else fractionAlong = 0.5;
	
	return true;
}

bool ParseProjectionPrePost(const ImportLogger &log, const pugi::xml_node &eProj,
	const Network &net, Network::Projection &proj
){
	
	auto prePopName = eProj.attribute("presynapticPopulation").value();
	auto postPopName = eProj.attribute("postsynapticPopulation").value();
	proj.presynapticPopulation = net.populations.get_id(prePopName);
	proj.postsynapticPopulation = net.populations.get_id(postPopName);
	if(proj.presynapticPopulation < 0){
		// extension!
		proj.presynapticPopulation = net.event_readers.get_id(prePopName);
		if(proj.presynapticPopulation < 0){
			log.error(eProj, "presynaptic population %s not found", prePopName);
			return false;
		}
		else proj.presyn_type = Network::Projection::PresynType::EVENT_SERIES;
	}
	else{
		proj.presyn_type = Network::Projection::PresynType::CELL; // it is just a cell population
	}
	if(proj.postsynapticPopulation < 0){
		log.error(eProj, "postsynaptic population %s not found", postPopName);
		return false;
	}
	
	return true;
}
// validate attached segment id's
bool ValidateAttachedSegmentId(
	const ImportLogger &log, const pugi::xml_node &eConn,
	const Morphology *morph, const char *sSegAttrName, Int cellId,
	Int &seg_seq
){
	
	Int segId = 0;
	if(!OptionalNonnegativeAttribute(log, eConn, sSegAttrName, segId)) return false;

	
	if( morph ){
		seg_seq = morph->segments.getSequential(segId); //convert to sequential form
		if(seg_seq < 0){
			log.error(eConn, "%s %ld not present in cell %ld", sSegAttrName, segId, cellId);
			return false;
		}
		return true; // seg_seq is known
	}
	else{
		// artificial cell
		if( segId != 0 ){
			log.error(eConn, "%s %ld not present in artificial cell %ld; only segment id 0 exists", sSegAttrName, segId, cellId);
			return false;
		}
		seg_seq = 0; // artificial cell only has "fictitious" segment 0
		return true;
	}
	// both cases covered
}
bool ParseConnectionPre_EventSeries(const ImportLogger &log, const pugi::xml_node &eConn,
	const Network::EventSetReader &evr,
	const char *cell_attr, const char *seg_attr,
	Int &cell_seq, Int &seg_seq, Real &frac
){
	// validate cell instance id's
	const char *sCell = RequiredAttribute(log, eConn, cell_attr);
	if(!sCell) return false;
	if(!( StrToL(sCell, cell_seq) && 0 <= cell_seq && cell_seq < evr.instances )){
		log.error(eConn, "%s must be an integer between 0 and %ld", cell_attr, (long) evr.instances);
		return false;
	}
	
	const char *sPort = eConn.attribute(seg_attr).value();
	if(*sPort){
		seg_seq = evr.ports.get_id(sPort);
		if(seg_seq < 0){
			log.error(eConn, "%s not found in event reader", sPort);
			return false;
		}
	}
	else{
		if(evr.ports.size() == 1){
			seg_seq = 0;
		}
		else if( (seg_seq = evr.ports.get_id("spike")) >= 0 ){
			// done in "if"
		}
		else{
			log.error(eConn, "attribute %s is needed, since the event reader has multiple ports but none of them is 'spike'", cell_attr, (long) evr.instances);
			return false;
		}
	}
	frac = 0.5; // to make the fake loca valid
	return true;
}
bool ParseConnectionPreOrPost_Cell(const ImportLogger &log, const pugi::xml_node &eConn,
	const Network::Population &pop, const Morphology *morphOrNot,
	const char *cell_attr, const char *seg_attr, const char *frac_attr,
	Int &cell_seq, Int &seg_seq, Real &frac
){
	// validate cell instance id's
	auto *sCell = eConn.attribute(cell_attr).value();
	if(!( *sCell )){
		log.error(eConn, "connection must have %s", cell_attr);
		return false;
	}
	
	// TODO shouldn't this be just an integer?
	Int cellId;
	if( !ParseSynapseCellRef(sCell, cellId) ){
		log.error(eConn, "invalid path \"%s\" for %s", sCell, cell_attr);
		return false;
	}
	cell_seq = pop.instances.getSequential(cellId); //convert to sequential form
	if(cell_seq < 0){
		log.error(eConn, "%s %ld not present in population", cell_attr, cellId);
		return false;
	}
	
	if( !ValidateAttachedSegmentId(log, eConn, morphOrNot, seg_attr, cellId, seg_seq) ) return false;
	
	// validate fraction along
	if(eConn.attribute(frac_attr)){
		if( !ParseQuantity<Dimensionless>(log, eConn, frac_attr, frac) ) return false;
		if(!( 0 <= frac && frac <= 1.0 )){
			log.error(eConn, "%s not between 0 and 1", frac_attr);
			return false;
		}
	}
	else frac = 0.5;
	
	return true;
}

bool ParseLoggerBase(const ImportLogger &log, const pugi::xml_node &eLogger, Simulation::LoggerBase &logger, bool is_classic_type){
	
	if(is_classic_type){
		const char *sFilename = RequiredAttribute(log, eLogger, "fileName");
		if(!sFilename) return false;
		logger.source_url = sFilename;
		// LATER enforce not using funny things or url's...
	}
	else{
		if(!GetSchemeOrRelpathFromUrl(log, eLogger, "href", logger.source_url)) return false;
		// NB: validate url on backend, until the url format is standardized LATER (or even better, leave it up to the backend?)
	}
	
	const char *sFormat = eLogger.attribute("format").value();
	logger.data_format = sFormat; // NB validate on backend, until the url format is standardized LATER (or even better, leave it up to the backend?)
	
	return true;
}

//The evolving result data of an import in progress
struct ImportState{
	
	// the model itself containing the following -- it is needed for path resolution while parsing the simulation parameters! (could break Simulation out of the model if needed, LATER)
	Model &model; // aliases with all the following, to avoiding putting "model." all over
	
	// hash table for re-process-able XML entities such as detached biophysics;
	// they make no sense without an associated Morphology
	CollectionWithNames<pugi::xml_node> standalone_biophysics;
	
	//references to Model
	DimensionSet &dimensions;
	CollectionWithNames<ComponentType> &component_types;
	CollectionWithNames<ComponentInstance> &component_instances;
	
	CollectionWithNames<Morphology> &morphologies;
	std::vector<BiophysicalProperties> &biophysics;
	
	CollectionWithNames<IonSpecies> &ion_species;
	CollectionWithNames<ConcentrationModel> &conc_models;
	
	CollectionWithNames<IonChannel> &ion_channels;
	
	CollectionWithNames<CellType> &cell_types;
	
	CollectionWithNames<SynapticComponent> &synaptic_components;
	
	CollectionWithNames<InputSource> &input_sources;
	
	CollectionWithNames<Network> &networks;
	
	CollectionWithNames<Simulation> &simulations;
	
	Int &target_simulation;
	
	// helpers for interface matching
	
	// TODO cache type compatibility checks with a hash table
	bool CheckSynapticComponentWithCellTypes(
		const ImportLogger &log, const pugi::xml_node &eConn,
		const SynapticComponent &syncomp, const char *syncomp_name, 
		const CellType &bound_cell, const char *bound_cell_name,
		const CellType &peer_cell, const char *peer_cell_name
	) const {
		
		// In the NeuroML model, a synapse consists of three or four parts:
		// 
		// On uni-directional, spike-based synapses:
		// Pre cell component -> Post synapse component <-> Post cell component
		//
		// On more general synapses:
		// Pre synapse component-?-Post synapse component
		//           |          \ /        |
		//           |           X         |
		//           |          / \        |
		//    Pre cell component   Post cell component
		//
		// For convenience, let's say cell and synapse on the same side are 'coupled', and a pre-part is a 'peer' to a post- part (and vice versa).
		// 
		// Though in theory synaptic component are connected to each other, all the information they need
		// 	is located in the underlying cell components of their peers; so they effectively tap into the cells themselves to receive the information.
		// These interactions may get more formally defined, in a future NeuroML/LEMS version (e.g. for proper STDP synapses).
		// 
		// Some cell models are physically dimensional, others are pure numerical (though time still is dimensional!)
		// General rules:
		//   Cell components are either compartments of physical cells, or point neurons entirely described in LEMS.
		//   If a synaptic component exposes  current, the coupled cell component must receive it, in the appropriate dimension.
		//   If a synaptic component requires voltage, the coupled cell component must expose  it, in the appropriate dimension.
		//   If a synaptic component receives  spikes, the peer    cell component must emit them.
		//   If a synaptic component requires   Vpeer, the peer    cell component must expose voltage, in the appropriate dimension.
		
		// the four general rules
		Dimension syn_current_dimension;
		Dimension bound_cell_current_dimension;
		Dimension syn_voltage_dimension;
		Dimension bound_cell_voltage_dimension;
		Dimension syn_vpeer_dimension;
		Dimension peer_cell_voltage_dimension;
		
		// bound cell compatible wrt current
		if( syncomp.GetCurrentOutputAndDimension( component_types, syn_current_dimension ) ){
			if( !bound_cell.GetCurrentInputAndDimension( component_types, bound_cell_current_dimension ) ){
				log.error( eConn, "synapse %s exposes current but cell %s does not receive current", syncomp_name, bound_cell_name);
				return false;
			}
			if( syn_current_dimension != bound_cell_current_dimension ){
				log.error( eConn, "synapse %s exposes current as %s but cell %s receives current as %s", syncomp_name, dimensions.Stringify(syn_current_dimension).c_str(), bound_cell_name, dimensions.Stringify(bound_cell_current_dimension).c_str() );
				return false;
			}
		}
		// bound cell compatible wrt voltage
		if( syncomp.GetVoltageInputAndDimension( component_types, syn_voltage_dimension ) ){
			if( !bound_cell.GetVoltageExposureAndDimension( component_types, bound_cell_voltage_dimension ) ){
				log.error( eConn, "synapse %s requires voltage but cell %s does not expose voltage", syncomp_name, bound_cell_name);
				return false;
			}
			if( syn_voltage_dimension != bound_cell_voltage_dimension ){
				log.error( eConn, "synapse %s requires voltage as %s but cell %s exposes voltage as %s", syncomp_name, dimensions.Stringify(syn_voltage_dimension).c_str(), bound_cell_name, dimensions.Stringify(bound_cell_voltage_dimension).c_str() );
				return false;
			}
		}
		// peer cell compatible wrt vpeer
		if( syncomp.GetVpeerInputAndDimension( component_types, syn_vpeer_dimension ) ){
			if( !peer_cell.GetVoltageExposureAndDimension( component_types, peer_cell_voltage_dimension ) ){
				log.error( eConn, "synapse %s exposes current but cell %s does not receive current", syncomp_name, bound_cell_name);
				return false;
			}
			if( syn_vpeer_dimension != peer_cell_voltage_dimension ){
				log.error( eConn, "synapse %s requires Vpeer as %s but cell %s exposes voltage as %s", syncomp_name, dimensions.Stringify(syn_vpeer_dimension).c_str(), bound_cell_name, dimensions.Stringify(peer_cell_voltage_dimension).c_str() );
				return false;
			}
		}
		// peer cell compatible wrt spikes
		if( syncomp.HasSpikeIn( component_types ) ){
			if( !peer_cell.HasSpikeOut( component_types, input_sources ) ){
				log.error( eConn, "synapse %s receives spikes but cell %s does not emit spikes", syncomp_name, bound_cell_name);
				return false;
			}
		}
		
		return true;
	};
	
	bool CheckInputComponentWithCellType(
		const ImportLogger &log, const pugi::xml_node &eInputInstance,
		const InputSource &input, const char *input_name, 
		Int cell_type_seq
	) const {
		
		const CellType &cello = cell_types.get(cell_type_seq);
		const char *cello_name = cell_types.getName(cell_type_seq);
		
		// Inputs are similar to synaptic components, but they're tied to one specific cell compartment, not two.
		//
		// General rules:
		//   Cell components are either compartments of physical cells, or point neurons entirely described in LEMS.
		//   If an input component exposes  current, the coupled cell component must receive it, in the appropriate dimension.
		//   If an input component requires voltage, the coupled cell component must expose  it, in the appropriate dimension.
		//   If an input component emits     spikes, the peer    cell component must receive them.
		
		// the three general rules
		Dimension input_current_dimension;
		Dimension cello_current_dimension;
		Dimension input_voltage_dimension;
		Dimension cello_voltage_dimension;
		
		// cell compatible wrt current
		if( input.GetCurrentOutputAndDimension( component_types, synaptic_components, input_current_dimension ) ){
			if( !cello.GetCurrentInputAndDimension( component_types, cello_current_dimension ) ){
				log.error( eInputInstance, "input source %s exposes current but cell %s does not receive current", input_name, cello_name);
				return false;
			}
			if( input_current_dimension != cello_current_dimension ){
				log.error( eInputInstance, "input source %s exposes current as %s but cell %s receives current as %s", input_name, dimensions.Stringify(input_current_dimension).c_str(), cello_name, dimensions.Stringify(cello_current_dimension).c_str() );
				return false;
			}
		}
		// cell compatible wrt voltage
		if( input.GetVoltageInputAndDimension( component_types, synaptic_components, input_voltage_dimension ) ){
			if( !cello.GetVoltageExposureAndDimension( component_types, cello_voltage_dimension ) ){
				log.error( eInputInstance, "input source %s requires voltage but cell %s does not expose voltage", input_name, cello_name);
				return false;
			}
			if( input_voltage_dimension != cello_voltage_dimension ){
				log.error( eInputInstance, "input source %s requires voltage as %s but cell %s exposes voltage as %s", input_name, dimensions.Stringify(input_voltage_dimension).c_str(), cello_name, dimensions.Stringify(cello_voltage_dimension).c_str() );
				return false;
			}
		}
		// peer cell compatible wrt spikes
		if( input.HasSpikeOut( component_types ) ){
			// remove check, because some core components do emit spikes as a side effect; reconsider LATER 
			if( 0 && !cello.HasSpikeIn( component_types ) ){
				log.error( eInputInstance, "input source %s emits spikes but cell %s does not receive spikes", input_name, cello_name);
				return false;
			}
		}
		
		return true;
	};
	
	// helpers for LEMS
	
	template<typename... Ts> struct make_void { typedef void type;};
	template<typename... Ts> using void_t = typename make_void<Ts...>::type;
	
	template <typename T, typename = void>
	struct is_iterable : std::false_type {};
	template <typename T>
	struct is_iterable<T, void_t<decltype(std::begin(std::declval<T>()))> > : std::true_type {};
	
	// LEMSify means: to make a ComponentInstance using parameters for built-in objects, following an already loaded ComponentType.
	// It is useful for giving the fully general LEMS treatment to instances of core NeuroML components, when there is no time to specialise each core NeuroML case (like all types of artifical cells)
	bool TryLemsifyParameter (
		const ImportLogger &log, const pugi::xml_node &eThing,
		const ComponentType &comptype, const char *prop_name, Real prop_value, ComponentInstance &compinst
	){
		Int prop_seq = comptype.properties.get_id(prop_name);
		if( prop_seq < 0 ){
			log.error(eThing, "internal error: lemsified property %s missing", prop_name);
			return false;
		}
		
		compinst.parms.push_back( { prop_seq, prop_value } );
		return true;
	};
	template< typename Functor, typename std::enable_if< !is_iterable<Functor>::value, int >::type = 7 >
	bool TryLemsifyComponent (
		const ImportLogger &log, const pugi::xml_node &eThing,
		const char *type_name, Functor DoParameterCustomizations, ComponentInstance &compinst
	){
		Int &comp_seq = compinst.id_seq;
		comp_seq = component_types.get_id(type_name);
		if( comp_seq < 0 ){
			log.error(eThing, "internal error: missing lemsified functor type %s", type_name);
			return false;
		}
		const auto &comp = component_types.get(comp_seq);
		return DoParameterCustomizations( log, eThing, comp, compinst );
	}
	bool TryLemsifyComponent (
		const ImportLogger &log, const pugi::xml_node &eThing,
		const char *type_name, ComponentInstance &compinst
	){
		return TryLemsifyComponent( log, eThing, type_name, []( auto &log, auto &eThing, auto &comp, auto &compinst ){ return true; }, compinst );
	}
	struct ParmEntry{
		const char *name;
		Real value; // in engine units
	};
	
	// lemsify component using each parameter
	template< typename Container >
	bool TryLemsifyComponent_WithParmContainer( const ImportLogger &log, const pugi::xml_node &eThing, const char *type_name, const Container &parms, ComponentInstance &compinst ){
		return TryLemsifyComponent(
			log, eThing, type_name,
			[ &parms, this ]( auto &log, auto &eThing, auto &comp, auto &compinst ){
				for( const auto &parm : parms ){
					if( !TryLemsifyParameter( log, eThing, comp, parm.name, parm.value, compinst ) ) return false;
				}
				return true;
			},
			compinst
		);
	}
	
	// for explicitly formed containers to pass
	template< size_t ParmCount >	
	bool TryLemsifyComponent ( const ImportLogger &log, const pugi::xml_node &eThing, const char *type_name, const  ParmEntry (&parms)[ParmCount], ComponentInstance &compinst ){
		return TryLemsifyComponent_WithParmContainer( log, eThing, type_name, parms, compinst );
	}
	
	// for brace-initialized, implicitly formed arrays to pass	
	template< typename Container, typename std::enable_if< is_iterable<Container>::value, int >::type = 7 >
	bool TryLemsifyComponent ( const ImportLogger &log, const pugi::xml_node &eThing, const char *type_name, const Container &parms, ComponentInstance &compinst ){
		return TryLemsifyComponent_WithParmContainer( log, eThing, type_name, parms, compinst );
	}
	
	// parse the tags one by one
	bool ParseStandaloneMorphology(const ImportLogger &log, const pugi::xml_node &eMorph){
		
		Morphology new_morph;
		
		auto name = RequiredNmlId(log, eMorph);
		if(!name) return false;
		
		if( !ParseMorphology(log, eMorph, new_morph)) return false;
		
		// add new standalone morphology, yay!
		morphologies.add(new_morph, name);
		return true;
	}
	
	bool ParseIonChannel(const ImportLogger &log, const pugi::xml_node &eChannel){
		
		IonChannel new_channel;
		
		auto name = RequiredNmlId(log, eChannel);
		if(!name) return false;
		//printf("Parsing ion channel %s\n",name);
		if( !::ParseIonChannel(log, eChannel, component_types, dimensions, ion_species, new_channel)) return false;
		
		// add new standalone channel, yay!
		ion_channels.add(new_channel, name);
		return true;
	}
	
	bool ParseStandaloneBiophysics(const ImportLogger &log, const pugi::xml_node &eBioph){
		
		auto name = RequiredNmlId(log, eBioph);
		if(!name) return false;
		
		// add new "standalone" biophysics, yay!
		standalone_biophysics.add(eBioph, name);
		return true;
	}
	
	bool ParsePhysicalCell(const ImportLogger &log, const pugi::xml_node &eCell, bool two_ca_pools = false){
		
		CellType cell_type;
		cell_type.type = CellType::PHYSICAL;
		
		auto name = RequiredNmlId(log, eCell);
		if(!name) return false;
		// printf("Parsing physical cell %s...\n", name);
		//parse cell properties
		//Required: id, also morph and bioph, internally or externally defined
		
		PhysicalCell &cell = cell_type.physical;
		
		//First get Morphology, as it is required for Biophysical properties to make sense
		//but morphology might as well be defined externally, so a symbolic form should exist too 
		Int morph_id = -1;
		if(eCell.child("morphology")){
			Morphology morph;
			if( !ParseMorphology(log, eCell.child("morphology"), morph)) return false;
			
			//add new Morphology to Morphologies array, yay! 
			morph_id = morphologies.add(morph);
		}
		else if(eCell.attribute("morphology")){
			//check externally defined Morphology
			const char *morph_name = eCell.attribute("morphology").value();
			morph_id = morphologies.get_id(morph_name);
			
			if(morph_id < 0){
				log.error(eCell," morphology attribute \"%s\" not found", morph_name);
				return false;
			}
		}
		else{
			log.error(eCell, "%s %s lacking morphology", eCell.name(), name );
			return false;
		}
		
		cell.morphology = morph_id;
		
		Int bioph_id = -1;
		const char * biophElementName = "biophysicalProperties";
		if(two_ca_pools) biophElementName = "biophysicalProperties2CaPools";
		
		pugi::xml_node eBioph;
		
		if(eCell.child(biophElementName)){
			//parse inline definition
			eBioph = eCell.child(biophElementName);
		}
		else if(eCell.attribute(biophElementName) && *eCell.attribute(biophElementName).value()){
			// check externally defined BiophysicalProperties, and realize them upon ths cell's morphology
			// clone the same Biophysical properties once for each cell type, this is probably not a problem
			auto biophStr = eCell.attribute(biophElementName).value();
			if(!standalone_biophysics.has(biophStr)){
				log.error(eCell, "biophysics type \"%s\" not found", biophStr );
				return false;
			}
			
			eBioph = standalone_biophysics.get(biophStr);
		}
		else{
			log.error(eCell, "%s %s lacking %s", eCell.name(), name, biophElementName );
			return false;
		}
		
		BiophysicalProperties bioph;
		
		if( !ParseBiophysicalProperties(log, eBioph, morphologies.get(morph_id), conc_models, ion_channels, component_types, two_ca_pools, ion_species, bioph)) return false;
		
		// ignore Standalone child elements
		
		//add new Biophysics to Biophysics array, yay! 
		biophysics.push_back(bioph);
		bioph_id = biophysics.size()-1;
		
		cell.biophysicalProperties = bioph_id;
		
		//add cell to cell types! yay!!
		cell_types.add(cell_type, name);
		return true;
	}
	
	bool ParseArtificialCell(const ImportLogger &log, const pugi::xml_node &eCell){
		
		CellType cell_type;
		cell_type.type = CellType::ARTIFICIAL;
		ArtificialCell &cell = cell_type.artificial;
		
		auto name = RequiredNmlId(log, eCell);
		if(!name) return false;
		
		const char *sType = TagNameOrType(eCell);
		
		// will be checked for unique name after parsing, try this approach anticipating conditional namespaces
		
		const static NameMap<ArtificialCell::Type> fakecell_types = {
			{"iafCell"                  , ArtificialCell::IAF                           },
			{"iafRefCell"               , ArtificialCell::IAF_REF                       },
			{"iafTauCell"               , ArtificialCell::IAF_TAU                       },
			{"iafTauRefCell"            , ArtificialCell::IAF_TAU_REF                   },
			{"izhikevichCell"           , ArtificialCell::IZH                           },
			{"izhikevich2007Cell"       , ArtificialCell::IZH_2007                      },
			{"adExIaFCell"              , ArtificialCell::ADEX                          },
			{"fitzHughNagumoCell"       , ArtificialCell::FN                            },
			{"fitzHughNagumo1969Cell"   , ArtificialCell::FN_1969                       },
			{"pinskyRinzelCA3Cell"      , ArtificialCell::PINSKY_RINZEL_CA3             },
			{"IF_curr_alpha"            , ArtificialCell::PYNN_IF_CURR_ALPHA            },
			{"IF_curr_exp"              , ArtificialCell::PYNN_IF_CURR_EXP              },
			{"IF_cond_alpha"            , ArtificialCell::PYNN_IF_COND_ALPHA            },
			{"IF_cond_exp"              , ArtificialCell::PYNN_IF_COND_EXP              },
			{"EIF_cond_exp_isfa_ista"   , ArtificialCell::PYNN_EIF_COND_EXP_ISFA_ISTA   },
			{"EIF_cond_alpha_isfa_ista" , ArtificialCell::PYNN_EIF_COND_ALPHA_ISFA_ISTA },
			{"HH_cond_exp"              , ArtificialCell::PYNN_HH_COND_EXP              },
			// {"compositeInput"           , ArtificialCell::COMPOSITE                     },
		};
		
		// first check if it's one of the core artificial cell types
		auto celltype_it = fakecell_types.find( sType );
		if(celltype_it != fakecell_types.end()){
			
			cell.type = celltype_it->second;
			
			auto ParseBasePynn = [ this ]( auto &log, const auto &eCell, ArtificialCell &cell ){
				if( !ParseQuantity<Dimensionless>(log, eCell, "cm"        , cell.cm        ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "i_offset"  , cell.i_offset  ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "tau_syn_E" , cell.tau_syn_E ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "tau_syn_I" , cell.tau_syn_I ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "v_init"    , cell.v_init    ) ) return false;
				return true;
			};
			auto ParseBasePynnIaf = [ ParseBasePynn, this ]( auto &log, const auto &eCell, ArtificialCell &cell ){
				if( !ParseBasePynn(log, eCell, cell) ) return false;
				
				if( !ParseQuantity<Dimensionless>(log, eCell, "tau_m"      , cell.tau_m   ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "tau_refrac" , cell.refract ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "v_reset"    , cell.reset   ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "v_rest"     , cell.v_rest  ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "v_thresh"   , cell.thresh  ) ) return false;
				
				return true;
			};
			auto ParseBasePynnIafCond = [ ParseBasePynnIaf, this ]( auto &log, const auto &eCell, ArtificialCell &cell ){
				if( !ParseBasePynnIaf(log, eCell, cell) ) return false;
				
				if( !ParseQuantity<Dimensionless>(log, eCell, "e_rev_E"    , cell.e_rev_E ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "e_rev_I"    , cell.e_rev_I ) ) return false;
				
				return true;
			};
			
			auto LemsifyParms_BasePynn = [ this ]( const auto &cell ){
				std::vector< ParmEntry > ret;
				ret.push_back( { "cm"        , cell.cm        } );
				ret.push_back( { "i_offset"  , cell.i_offset  } );
				ret.push_back( { "tau_syn_E" , cell.tau_syn_E } );
				ret.push_back( { "tau_syn_I" , cell.tau_syn_I } );
				ret.push_back( { "v_init"    , cell.v_init    } );
				return ret;
			};
			auto LemsifyParms_BasePynnIaf = [ LemsifyParms_BasePynn, this ]( const auto &cell ){
				std::vector< ParmEntry > ret = LemsifyParms_BasePynn(cell);
				ret.push_back( { "tau_m"      , cell.tau_m   } );
				ret.push_back( { "tau_refrac" , cell.refract } );
				ret.push_back( { "v_reset"    , cell.reset   } );
				ret.push_back( { "v_rest"     , cell.v_rest  } );
				ret.push_back( { "v_thresh"   , cell.thresh  } );
				return ret;
			};
			auto LemsifyParms_BasePynnIafCond = [ LemsifyParms_BasePynnIaf, this ]( const auto &cell ){
				std::vector< ParmEntry > ret = LemsifyParms_BasePynnIaf(cell);
				ret.push_back( { "e_rev_E"    , cell.e_rev_E } );
				ret.push_back( { "e_rev_I"    , cell.e_rev_I } );
				return ret;
			};
			
			if( cell.type == ArtificialCell::IAF || cell.type == ArtificialCell::IAF_REF ){
				if( !ParseQuantity<Capacitance>(log, eCell, "C", cell.C ) ) return false;
				if( !ParseQuantity<Conductance>(log, eCell, "leakConductance", cell.leakConductance ) ) return false;
				if( !ParseQuantity<Voltage    >(log, eCell, "leakReversal"   , cell.leakReversal    ) ) return false;
				if( !ParseQuantity<Voltage>(log, eCell, "thresh", cell.thresh ) ) return false;
				if( !ParseQuantity<Voltage>(log, eCell, "reset" , cell.reset  ) ) return false;
				if( cell.type == ArtificialCell::IAF_REF ){
					if( !ParseQuantity<Time>(log, eCell, "refract" , cell.refract  ) ) return false;
				}
			}
			else if( cell.type == ArtificialCell::IAF_TAU || cell.type == ArtificialCell::IAF_TAU_REF ){
				if( !ParseQuantity<Time>(log, eCell, "tau", cell.tau ) ) return false;
				if( !ParseQuantity<Voltage    >(log, eCell, "leakReversal"   , cell.leakReversal    ) ) return false;
				if( !ParseQuantity<Voltage>(log, eCell, "thresh", cell.thresh ) ) return false;
				if( !ParseQuantity<Voltage>(log, eCell, "reset" , cell.reset  ) ) return false;
				if( cell.type == ArtificialCell::IAF_TAU_REF ){
					if( !ParseQuantity<Time>(log, eCell, "refract" , cell.refract  ) ) return false;
				}
			}
			else if( cell.type == ArtificialCell::IZH ){
				if( !ParseQuantity<Voltage>(log, eCell, "v0"     , cell.v0  ) ) return false;
				if( !ParseQuantity<Voltage>(log, eCell, "thresh" , cell.thresh ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "a", cell.a ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "b", cell.b ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "c", cell.c ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "d", cell.d ) ) return false;
			}
			else if( cell.type == ArtificialCell::IZH_2007 ){
				if( !ParseQuantity<Capacitance>(log, eCell, "C", cell.C ) ) return false;
				
				if( !ParseQuantity<Voltage>(log, eCell, "v0"    , cell.v0     ) ) return false;
				if( !ParseQuantity<Voltage>(log, eCell, "vr"    , cell.reset  ) ) return false;
				if( !ParseQuantity<Voltage>(log, eCell, "vt"    , cell.thresh ) ) return false;
				if( !ParseQuantity<Voltage>(log, eCell, "vpeak" , cell.vpeak  ) ) return false;
				if( !ParseLemsQuantity(log, eCell, "k" , dimensions, dimensions.Get("conductance_per_voltage"), cell.k  ) ) return false;
				if( !ParseQuantity<Frequency  >(log, eCell, "a", cell.a ) ) return false;
				if( !ParseQuantity<Conductance>(log, eCell, "b", cell.b ) ) return false;
				if( !ParseQuantity<Voltage    >(log, eCell, "c", cell.c ) ) return false;
				if( !ParseQuantity<Current    >(log, eCell, "d", cell.d ) ) return false;
			}
			else if( cell.type == ArtificialCell::ADEX ){
				if( !ParseQuantity<Capacitance>(log, eCell, "C", cell.C ) ) return false;
				
				if( !ParseQuantity<Conductance>(log, eCell, "gL", cell.leakConductance ) ) return false;
				if( !ParseQuantity<Voltage    >(log, eCell, "EL", cell.leakReversal    ) ) return false;
				
				if( !ParseQuantity<Voltage>(log, eCell, "thresh", cell.thresh ) ) return false;
				if( !ParseQuantity<Voltage>(log, eCell, "reset" , cell.reset  ) ) return false;
				if( !ParseQuantity<Voltage>(log, eCell, "VT"    , cell.vt     ) ) return false;
				if( !ParseQuantity<Voltage>(log, eCell, "delT"  , cell.delt   ) ) return false;
				
				if( !ParseQuantity<Time>(log, eCell, "tauw"    , cell.tau      ) ) return false;
				if( !ParseQuantity<Time>(log, eCell, "refract" , cell.refract  ) ) return false;
				if( !ParseQuantity<Conductance>(log, eCell, "a", cell.a ) ) return false;
				if( !ParseQuantity<Current    >(log, eCell, "b", cell.b ) ) return false;
			}
			else if( cell.type == ArtificialCell::FN ){
				if( !ParseQuantity<Dimensionless>(log, eCell, "I", cell.I ) ) return false;
			}
			else if( cell.type == ArtificialCell::FN_1969 ){
				if( !ParseQuantity<Dimensionless>(log, eCell, "I", cell.I ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "a", cell.a ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "b", cell.b ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "phi", cell.phi ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "V0", cell.v0 ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "W0", cell.w0 ) ) return false;
			}
			else if( cell.type == ArtificialCell::PINSKY_RINZEL_CA3 ){
				if( !ParseLemsQuantity(log, eCell, "iSoma" , dimensions, dimensions.Get("currentDensity"), cell.iSoma  ) ) return false;
				if( !ParseLemsQuantity(log, eCell, "iDend" , dimensions, dimensions.Get("currentDensity"), cell.iDend  ) ) return false;
				if( !ParseQuantity<Conductivity  >(log, eCell, "gNmda"  , cell.gNmda  ) ) return false;
				if( !ParseQuantity<Conductivity  >(log, eCell, "gAmpa"  , cell.gAmpa  ) ) return false;
				if( !ParseQuantity<Conductivity  >(log, eCell, "gc"     , cell.gc     ) ) return false;
				if( !ParseQuantity<Conductivity  >(log, eCell, "gLs"    , cell.gLs    ) ) return false;
				if( !ParseQuantity<Conductivity  >(log, eCell, "gLd"    , cell.gLd    ) ) return false;
				if( !ParseQuantity<Conductivity  >(log, eCell, "gNa"    , cell.gNa    ) ) return false;
				if( !ParseQuantity<Conductivity  >(log, eCell, "gKdr"   , cell.gKdr   ) ) return false;
				if( !ParseQuantity<Conductivity  >(log, eCell, "gCa"    , cell.gCa    ) ) return false;
				if( !ParseQuantity<Conductivity  >(log, eCell, "gKahp"  , cell.gKahp  ) ) return false;
				if( !ParseQuantity<Conductivity  >(log, eCell, "gKC"    , cell.gKC    ) ) return false;
				if( !ParseQuantity<Voltage       >(log, eCell, "eNa"    , cell.eNa    ) ) return false;
				if( !ParseQuantity<Voltage       >(log, eCell, "eCa"    , cell.eCa    ) ) return false;
				if( !ParseQuantity<Voltage       >(log, eCell, "eK"     , cell.eK     ) ) return false;
				if( !ParseQuantity<Voltage       >(log, eCell, "eL"     , cell.eL     ) ) return false;
				if( !ParseQuantity<Dimensionless >(log, eCell, "qd0"    , cell.qd0    ) ) return false;
				if( !ParseQuantity<Dimensionless >(log, eCell, "pp"     , cell.pp     ) ) return false;
				if( !ParseQuantity<Dimensionless >(log, eCell, "alphac" , cell.alphac ) ) return false;
				if( !ParseQuantity<Dimensionless >(log, eCell, "betac"  , cell.betac  ) ) return false;
				
				if( !ParseQuantity<SpecificCapacitance>(log, eCell, "cm"  , cell.cm  ) ) return false;
			}
			else if( cell.type == ArtificialCell::PYNN_IF_CURR_ALPHA || cell.type == ArtificialCell::PYNN_IF_CURR_EXP ){
				if( !ParseBasePynnIaf(log, eCell, cell ) ) return false;
			}
			else if( cell.type == ArtificialCell::PYNN_IF_COND_ALPHA || cell.type == ArtificialCell::PYNN_IF_COND_EXP ){
				if( !ParseBasePynnIafCond(log, eCell, cell ) ) return false;
			}
			else if( cell.type == ArtificialCell::PYNN_EIF_COND_EXP_ISFA_ISTA || cell.type == ArtificialCell::PYNN_EIF_COND_ALPHA_ISFA_ISTA ){
				if( !ParseBasePynnIafCond(log, eCell, cell ) ) return false;
				
				if( !ParseQuantity<Dimensionless>(log, eCell, "a"       , cell.a       ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "b"       , cell.b       ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "delta_T" , cell.delt    ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "tau_w"   , cell.tau_w   ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "v_spike" , cell.v_spike ) ) return false;
				
			}
			else if( cell.type == ArtificialCell::PYNN_HH_COND_EXP ){
				if( !ParseBasePynn(log, eCell, cell ) ) return false;
				
				if( !ParseQuantity<Dimensionless>(log, eCell, "v_offset"   , cell.v_offset   ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "e_rev_E"    , cell.e_rev_E    ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "e_rev_I"    , cell.e_rev_I    ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "e_rev_K"    , cell.e_rev_K    ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "e_rev_Na"   , cell.e_rev_Na   ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "e_rev_leak" , cell.e_rev_leak ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "g_leak"     , cell.g_leak     ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "gbar_K"     , cell.gbar_K     ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eCell, "gbar_Na"    , cell.gbar_Na    ) ) return false;
				
			}
			else{
				log.error(eCell, "internal error: unknown artificial cell type %s", sType );
				return false;
			}

			// along with the LEMS component, attach a mapping to the LEMS version
			// also initialize the id_seq to show it's missing
			cell.component.clear();
			
			if( cell.type == ArtificialCell::IAF ){
				if( !TryLemsifyComponent( log, eCell, sType, {
					{ "C" , cell.C  },
					{ "leakConductance", cell.leakConductance },
					{ "leakReversal"   , cell.leakReversal    },
					{ "thresh", cell.thresh },
					{ "reset" , cell.reset  },
				}, cell.component ) ) return false;
			}
			else if( cell.type == ArtificialCell::IAF_REF ){
				if( !TryLemsifyComponent( log, eCell, sType, {
					{ "C" , cell.C  },
					{ "leakConductance", cell.leakConductance },
					{ "leakReversal"   , cell.leakReversal    },
					{ "thresh", cell.thresh },
					{ "reset" , cell.reset  },
					{ "refract" , cell.refract  },
				}, cell.component ) ) return false;
			}
			else if( cell.type == ArtificialCell::IAF_TAU ){
				if( !TryLemsifyComponent( log, eCell, sType, {
					{ "tau" , cell.tau  },
					{ "leakReversal"   , cell.leakReversal    },
					{ "thresh", cell.thresh },
					{ "reset" , cell.reset  },
				}, cell.component ) ) return false;
			}
			else if( cell.type == ArtificialCell::IAF_TAU_REF ){
				if( !TryLemsifyComponent( log, eCell, sType, {
					{ "tau" , cell.tau  },
					{ "leakReversal"   , cell.leakReversal    },
					{ "thresh", cell.thresh },
					{ "reset" , cell.reset  },
					{ "refract" , cell.refract  },
				}, cell.component ) ) return false;
			}
			else if( cell.type == ArtificialCell::IZH ){
				if( !TryLemsifyComponent( log, eCell, sType, {
					{ "v0" , cell.v0  },
					{ "thresh", cell.thresh },
					{ "a" , cell.a  },
					{ "b" , cell.b  },
					{ "c" , cell.c  },
					{ "d" , cell.d  },
				}, cell.component ) ) return false;
			}
			else if( cell.type == ArtificialCell::IZH_2007 ){
				if( !TryLemsifyComponent( log, eCell, sType, {
					{ "C" , cell.C  },
					{ "v0" , cell.v0  },
					{ "vr", cell.reset },
					{ "vt", cell.thresh },
					{ "vpeak", cell.vpeak },
					{ "k" , cell.k  },
					{ "a" , cell.a  },
					{ "b" , cell.b  },
					{ "c" , cell.c  },
					{ "d" , cell.d  },
				}, cell.component ) ) return false;
			}
			else if( cell.type == ArtificialCell::IZH_2007 ){
				if( !TryLemsifyComponent( log, eCell, sType, {
					{ "C" , cell.C  },
					{ "v0" , cell.v0  },
					{ "vr", cell.reset },
					{ "vt", cell.thresh },
					{ "vpeak", cell.vpeak },
					{ "k" , cell.k  },
					{ "a" , cell.a  },
					{ "b" , cell.b  },
					{ "c" , cell.c  },
					{ "d" , cell.d  },
				}, cell.component ) ) return false;
			}
			else if( cell.type == ArtificialCell::ADEX ){
				if( !TryLemsifyComponent( log, eCell, sType, {
					{ "C" , cell.C },
					{ "gL" , cell.leakConductance },
					{ "EL", cell.leakReversal },
					{ "thresh", cell.thresh },
					{ "reset", cell.reset },
					{ "VT" , cell.vt },
					{ "delT" , cell.delt },
					{ "tauw" , cell.tau },
					{ "refract" , cell.refract },
					{ "a" , cell.a  },
					{ "b" , cell.b  },
				}, cell.component ) ) return false;
			}
			else if( cell.type == ArtificialCell::FN ){
				if( !TryLemsifyComponent( log, eCell, sType, {
					{ "I" , cell.I },
				}, cell.component ) ) return false;
			}
			else if( cell.type == ArtificialCell::FN_1969 ){
				if( !TryLemsifyComponent( log, eCell, sType, {
					{ "I" , cell.I },
					{ "a" , cell.a },
					{ "b" , cell.b },
					{ "phi" , cell.phi },
					{ "V0" , cell.v0 },
					{ "W0" , cell.w0 },
				}, cell.component ) ) return false;
			}
			else if( cell.type == ArtificialCell::PINSKY_RINZEL_CA3 ){
				if( !TryLemsifyComponent( log, eCell, sType, {
					{ "iSoma"  , cell.iSoma  },
					{ "iDend"  , cell.iDend  },
					{ "gNmda"  , cell.gNmda  },
					{ "gAmpa"  , cell.gAmpa  },
					{ "gc"     , cell.gc     },
					{ "gLs"    , cell.gLs    },
					{ "gLd"    , cell.gLd    },
					{ "gNa"    , cell.gNa    },
					{ "gKdr"   , cell.gKdr   },
					{ "gCa"    , cell.gCa    },
					{ "gKahp"  , cell.gKahp  },
					{ "gKC"    , cell.gKC    },
					{ "eNa"    , cell.eNa    },
					{ "eCa"    , cell.eCa    },
					{ "eK"     , cell.eK     },
					{ "eL"     , cell.eL     },
					{ "qd0"    , cell.qd0    },
					{ "pp"     , cell.pp     },
					{ "alphac" , cell.alphac },
					{ "betac"  , cell.betac  },
					{ "cm"     , cell.cm     },
				}, cell.component ) ) return false;
			}
			else if( cell.type == ArtificialCell::PYNN_IF_CURR_ALPHA || cell.type == ArtificialCell::PYNN_IF_CURR_EXP ){
				if( !TryLemsifyComponent( log, eCell, sType, LemsifyParms_BasePynnIaf(cell), cell.component ) ) return false;
			}
			else if( cell.type == ArtificialCell::PYNN_IF_COND_ALPHA || cell.type == ArtificialCell::PYNN_IF_COND_EXP ){
				if( !TryLemsifyComponent( log, eCell, sType, LemsifyParms_BasePynnIafCond(cell), cell.component ) ) return false;
			}
			else if( cell.type == ArtificialCell::PYNN_EIF_COND_EXP_ISFA_ISTA || cell.type == ArtificialCell::PYNN_EIF_COND_ALPHA_ISFA_ISTA ){
				auto ret = LemsifyParms_BasePynnIafCond(cell);
				ret.push_back( { "a"       , cell.a       } );
				ret.push_back( { "b"       , cell.b       } );
				ret.push_back( { "delta_T" , cell.delt    } );
				ret.push_back( { "tau_w"   , cell.tau_w   } );
				ret.push_back( { "v_spike" , cell.v_spike } );
				if( !TryLemsifyComponent( log, eCell, sType, ret, cell.component ) ) return false;
			}
			else if( cell.type == ArtificialCell::PYNN_HH_COND_EXP ){
				auto ret = LemsifyParms_BasePynn(cell);
				ret.push_back( { "v_offset"   , cell.v_offset    } );
				ret.push_back( { "e_rev_E"    , cell.e_rev_E     } );
				ret.push_back( { "e_rev_I"    , cell.e_rev_I     } );
				ret.push_back( { "e_rev_K"    , cell.e_rev_K     } );
				ret.push_back( { "e_rev_Na"   , cell.e_rev_Na    } );
				ret.push_back( { "e_rev_leak" , cell.e_rev_leak  } );
				ret.push_back( { "g_leak"     , cell.g_leak      } );
				ret.push_back( { "gbar_K"     , cell.gbar_K      } );
				ret.push_back( { "gbar_Na"    , cell.gbar_Na     } );
				if( !TryLemsifyComponent( log, eCell, sType, ret, cell.component ) ) return false;
			}
			else{
				log.error(eCell, "internal error: unknown lemsify artificial cell" );
				return false;
			}
			
		}
		else{
			// should be a LEMS component
			// perhaps check parm completion for the LEMSified core components too? nah, let the NANs leak instead (it would be an internal bug anyway)
			
			if( !ParseComponentInstanceArtificialCell(log, eCell, component_types, dimensions, sType, cell.component) ) return false;
			
			cell.type = ArtificialCell::COMPONENT;
			// log.error(eCell, "unknown cell type %s", eCell.name() );
			// return false;
		}
		
		// now check if fully-specified cell type is new to its namespace
		if( cell_types.has(name) ){
			log.error(eCell, "cell type %s already defined", name );
			return false;
		}
		
		//add cell to artificial cells! yay!
		cell_types.add(cell_type, name);
		
		return true;
	}
	
	bool ParseConcentrationModel(const ImportLogger &log, const pugi::xml_node &ePool){
		
		ConcentrationModel new_pool;
		
		auto name = RequiredNmlId(log, ePool);
		if(!name) return false;
		
		const char *ion_name = ePool.attribute("ion").value();
		if( !*ion_name ){
			log.error(ePool, "ion concentration model missing ion species attribute");
			return false;
		}
		//convert ion species name to id
		new_pool.ion_species = ion_species.idOrNew(ion_name);
		
		auto type = TagNameOrType(ePool);
		// add a case when "concnetrationModel" remains unspecified
		if(strcmp(type, "concentrationModel") == 0 && *ePool.attribute("type").value() ){
			log.error( ePool, "<%s> does not specify type value", type );
			return false;
		}
		else if(strcmp(type, "decayingPoolConcentrationModel") == 0){
			new_pool.type = ConcentrationModel::LEAKY;
			if( !ParseQuantity<Length>(log, ePool, "shellThickness", new_pool.shellThickness_or_rhoFactor) ) return false;
			// common for decayingPoolConcentrationModel and fixedFactorConcentrationModel
			if( !ParseQuantity<Concentration>(log, ePool, "restingConc", new_pool.restingConc) ) return false;
			if( !ParseQuantity<Time>(log, ePool, "decayConstant", new_pool.decayConstant) ) return false;
		}
		else if(strcmp(type, "fixedFactorConcentrationModel") == 0){
			new_pool.type = ConcentrationModel::FIXED_FACTOR;
			if( !ParseQuantity<RhoFactor>(log, ePool, "rho", new_pool.shellThickness_or_rhoFactor) ) return false;
			// common for decayingPoolConcentrationModel and fixedFactorConcentrationModel
			if( !ParseQuantity<Concentration>(log, ePool, "restingConc", new_pool.restingConc) ) return false;
			if( !ParseQuantity<Time>(log, ePool, "decayConstant", new_pool.decayConstant) ) return false;
		}
		else{
			// must be a LEMS type then (it had to be detected as a compatible component type anyway, for control flow to be redirected here)
			// TODO add better error logging to indicate this, leading to specified "component type"
			//log.error(ePool, "how did a concentration model of type %s get here?", type);
			new_pool.type = ConcentrationModel::COMPONENT;
			if( !ParseComponentInstanceConcentrationModel(log, ePool, component_types, dimensions, type, new_pool.component) ) return false;
		}
		// ignore Standalone children: notes, annotation, property
		
		//add pool to concentration models! yay!
		conc_models.add(new_pool, name);
		
		return true;
	}
	
	// helpers here
	
	bool ParseNetwork(const ImportLogger &log, const pugi::xml_node &eNet){
		
		Network net;
		
		auto name = RequiredNmlId(log, eNet);
		if(!name) return false;
		//TODO redefinition sanity check!
		
		bool has_temperature = false; // must specify, though!
		bool needs_temperature = true; // should whine only when a temperature-dependent thing is actually present, though
		auto sType = eNet.attribute("type").value();
		if(*sType){
			if(strcmp(sType,"network") == 0){
				//well, duh
			}
			else if(strcmp(sType,"networkWithTemperature") == 0){
				has_temperature = true;
			}
			else{
				log.error(eNet, "unknown network type %s", sType);
				return false;
			}
		}
		if( strcmp(eNet.name(), "networkWithTemperature") == 0){
			has_temperature = true;
		}
		
		if( needs_temperature ){
			if( has_temperature ){
				if( !ParseQuantity<Temperature>(log, eNet, "temperature", net.temperature) ) return false;
			}
			else{
				// whine
				log.warning(eNet, "Temperature not specified, set to NEURON default: 6.3 degC");
				net.temperature = 6.3 + 273.15; // Kelvin degrees
			}
		}
		// TODO order them by element name maybe?
		for (auto eNetEl: eNet.children()){
			// printf("%s\n", eNetEl.name());
			// perhaps annotation stuff?
			// NOTE wherever there is Standalone, notes, annotation and property tage may exist. Handle properties as they appear, for they must be documented to be used.
			if(
				strcmp( eNetEl.name(), "notes" ) == 0
				|| strcmp( eNetEl.name(), "annotation" ) == 0
			){
				continue;
			}
			if(strcmp(eNetEl.name(), "space") == 0){
				// this tag is used precisely nowhere outside the XSD
				// use LATER when it is actually used somewhere
			}
			else if(strcmp(eNetEl.name(), "region") == 0){
				// this tag is used precisely nowhere outside the XSD
				// use LATER when it is actually used somewhere
			}
			else if(strcmp(eNetEl.name(), "extracellularProperties") == 0){
				// const auto &eExtra = eNetEl;
				// what do the species spec semantics mean here? apply to entire cell, or what?
				log.error(eNetEl, "extracellularProperties in network specification not supported yet");
				return false;
			}
			else if(strcmp(eNetEl.name(), "population") == 0){
				const auto &eUniPop = eNetEl;
				// an uniform population, consisting of clones of a particular cell
				Network::Population pop;
				
				auto pop_name = RequiredNmlId(log, eUniPop);
				if(!pop_name) return false;
				
				auto celltype_name = RequiredAttribute(log, eUniPop, "component");
				if(!celltype_name) return false;
				pop.component_cell = cell_types.get_id(celltype_name);
				if(pop.component_cell < 0){
					
					// not a cell type, but it could be a spiking input source !
					Int input_source_seq = input_sources.get_id(celltype_name);
					// full search in LEMS namespace, perhaps LATER
					if( input_source_seq >= 0 ){
						const auto &input = input_sources.get(input_source_seq);
						
						// try a spike input as artificial cell
						if( !RealizeInputSourceAsArtificialCell( log, eUniPop, input, celltype_name ) )return false;
						pop.component_cell = cell_types.get_id(celltype_name);
						if(pop.component_cell < 0){
							log.error( eUniPop, "internal error: could not create cell type from input source somehow, despite validating" );
							return false;
						}
					}
					else{
						//try realizing it as a component
						
						Int compinst_seq = component_instances.get_id(celltype_name);
						if( compinst_seq < 0 ){
							log.error( eUniPop, "unknown cell type %s", celltype_name );
							return false;
						}
						
						const auto &compinst = component_instances.get(compinst_seq);
						
						if( !RealizeComponentAsArtificialCell( log, eUniPop, compinst, celltype_name ) ) return false;
						// log.error(eUniPop, "component cell type %s not found in cell types or input sources", celltype_name);
						pop.component_cell = cell_types.get_id(celltype_name);
						if(pop.component_cell < 0){
							log.error( eUniPop, "internal error: could not create cell type from component type somehow, despite validating" );
							return false;
						}
					}
				}
				
				// but are these instances named and placed, or not?
				// let's try auto_detecting
				
				// Size and instances may both be defined, redundantly
				// Let specific instances take priority, but also validate with 'size' attribute if available
				if(eUniPop.child("instance")){
					for (auto ePopEl: eUniPop.children()){
						if(strcmp(ePopEl.name(), "instance") == 0){
							const auto &eInstance = ePopEl;
							Network::Population::Instance instance;
							Int instance_id;
							if( !StrToL(eInstance.attribute("id").value(), instance_id) ){
								log.error(eInstance, "instance requires a numeric id");
								return false;
							}
							if(instance_id < 0){
								log.error(eInstance, "instance id is negative");
								return false;
							}
							
							const auto &eLocation = eInstance.child("location");
							if(!eLocation){
								log.error(eInstance, "instance requires a <location> tag");
								return false;
							}
							
							//read 3d point of instnatiation in dimensionless(implied microns), same as Morphology coordinates
							if(!(
								   StrToF(eLocation.attribute("x").value(), instance.x)
								&& StrToF(eLocation.attribute("y").value(), instance.y)
								&& StrToF(eLocation.attribute("z").value(), instance.z)
							)){
								log.error(eLocation,
									"instance %ld has an invalid distal point (%s, %s, %s) position", instance_id,
									eLocation.attribute("x").value(),
									eLocation.attribute("y").value(),
									eLocation.attribute("z").value()
								);
								return false;
							}
							
							pop.instances.add(instance, instance_id);
							
						}
						else if(strcmp(ePopEl.name(), "layout") == 0){
							// actually implemented in the JSON-based NeuroMLLite spec
							log.error(ePopEl, "layout not supported yet");
						}
						else{
							// unknown, skip
							// also ignore Standalone child elements
						}
					}
				}
				// if instances were already defined, validate size; else synthesize a collection of cells
				if(eUniPop.attribute("size")){
					Int size;
					if(!RequiredPositiveAttribute(log, eUniPop, "size", size)) return false;
					if( !pop.instances.contents.empty() ){
						if( (int)pop.instances.contents.size() != size){
							log.error(eUniPop, "Population size is %zd, while explicit instances were %ld", pop.instances.contents.size(), size);
							return false;
						}
					}
					else{
						for( Int i = 0; i < size; i++ ){
							Network::Population::Instance instance = {0.0, 0.0, 0.0};
							pop.instances.add(instance, i);
						}
					}
				}
				
				// add population to network! yay!
				net.populations.add(pop, pop_name);
			}
			else if(strcmp(eNetEl.name(), "cellSet") == 0){
				// const auto &eSet = eNetEl;
				// nobody has actually implemented this
				log.error(eNetEl, "cellSet not supported yet");
				return false;
			}
			else if(strcmp(eNetEl.name(), "synapticConnection") == 0){
				// const auto &eConn = eNetEl;
				// this is just for LEMS to work, but aimed toward artifical cells
				// exists in NeuroML examples, reason enough to support it
				
				
				log.error(eNetEl, "synapticConnection not supported yet");
				return false;
			}
			else if(
				   strcmp(eNetEl.name(), "projection") == 0
				|| strcmp(eNetEl.name(), "electricalProjection") == 0
				|| strcmp(eNetEl.name(), "continuousProjection") == 0
			){
				const auto &eProj = eNetEl;
				
				Network::Projection proj;
				
				auto name = RequiredNmlId(log, eProj);
				if(!name) return false;
				
				bool is_spiking = (strcmp(eNetEl.name(), "projection") == 0); //projection is supposed to be for spiking connections
				//electricalProjection is supposed to be for linear gap junctions
				//continuousProjection is supposed to be for any graded synapse, why not gap junctions too?
				
				if(!ParseProjectionPrePost(log, eProj, net, proj)) return false;
				
				Int default_synapse_type = -1;
				if(is_spiking){
					// may be present for all projections LATER
					auto synName = eProj.attribute("synapse").value();
					if( !*synName ){
						log.error(eProj, "projection requires synapse type", synName);
						return false;
					}
					default_synapse_type = synaptic_components.get_id(synName);
					if(default_synapse_type < 0){
						log.error(eProj, "projection synapse type %s not found", synName);
						return false;
					}
				}
				
				// get morphologies, they will be needed to validate cell segment - cell segment connections
				// and get some nullable properties, or rather std::optional?
				const Network::EventSetReader *presynaptic_evr = NULL;
				const Network::Population *presynaptic_population = NULL;
				const CellType *cell_type_pre = NULL;
				const char *cell_type_name_pre = NULL;
				const Morphology *pre_morph  = NULL;
				if( proj.presyn_type == Network::Projection::PresynType::EVENT_SERIES ){
					presynaptic_evr = &net.event_readers.get(proj.presynapticPopulation);
					cell_type_name_pre = net.event_readers.getName(proj.presynapticPopulation);
				}
				else if( proj.presyn_type == Network::Projection::PresynType::CELL ){
					presynaptic_population = &net.populations.get(proj.presynapticPopulation);
					cell_type_pre = &cell_types.get( presynaptic_population->component_cell );
					cell_type_name_pre = cell_types.getName( presynaptic_population->component_cell );
					if( cell_type_pre->type == CellType::PHYSICAL ) pre_morph  = &( morphologies.get( cell_type_pre->physical.morphology ) );
				}
				else{ assert(false); return false; }
				
				const Network::Population &postsynaptic_population = net.populations.get(proj.postsynapticPopulation);
				const CellType &cell_type_post = cell_types.get( postsynaptic_population.component_cell );
				const char *cell_type_name_post = cell_types.getName( postsynaptic_population.component_cell );
				const Morphology *post_morph = NULL;
				if( cell_type_post.type == CellType::PHYSICAL ) post_morph = &( morphologies.get( cell_type_post.physical.morphology ) );
				
				// and also determine properties of subtypes
				struct ConnectionSubType{
					Network::Projection::Connection::Type type;
					bool uses_weight;
					bool uses_delay;
				};
				const static NameMap<ConnectionSubType> connection_types = {
					{"connection"					, {Network::Projection::Connection::SPIKING		,false, false} },
					{"connectionWD"					, {Network::Projection::Connection::SPIKING		,true , true } },
					{"electricalConnection"			, {Network::Projection::Connection::ELECTRICAL	,false, false} },
					{"electricalConnectionInstance"	, {Network::Projection::Connection::ELECTRICAL	,false, false} },
					{"electricalConnectionInstanceW", {Network::Projection::Connection::ELECTRICAL	,true , false} },
					{"continuousConnection"			, {Network::Projection::Connection::CONTINUOUS	,false, false} },
					{"continuousConnectionInstance"	, {Network::Projection::Connection::CONTINUOUS	,false, false} },
					{"continuousConnectionInstanceW", {Network::Projection::Connection::CONTINUOUS	,true , false} },
				};
				
				for(auto eProjEl: eProj.children()){
					// first check if it's one of the connection types
					auto conntype_it = connection_types.find(eProjEl.name());
					if(conntype_it != connection_types.end()){
						//printf("%s %s\n", eNetEl.name(), eProjEl.name());
						const auto &eConn = eProjEl;
						Network::Projection::Connection conn;
						
						conn.type = conntype_it->second.type;
						bool uses_weight = conntype_it->second.uses_weight;
						bool uses_delay = conntype_it->second.uses_delay;
						bool uses_old_format = ( conn.type == Network::Projection::Connection::SPIKING );
						
						if( proj.presyn_type == Network::Projection::PresynType::EVENT_SERIES && conn.type != Network::Projection::Connection::SPIKING ){
							log.error(eConn, "projection from event series must have one-way synapses only");
							return false;
						}
						
						Int id;
						if(!( StrToL(eConn.attribute("id").value(), id) )){
							log.error(eConn, "connection requires a non-negative id");
							return false;
						}
						if(id < 0){
							log.error(eConn, "connection id is negative");
							return false;
						}
						if(proj.connections.hasId(id)){
							log.error(eConn, "connection %ld already defined", id);
							return false;
						}
						
						// get synapse type from XML element and check corresponding rules
						
						// Projection rules:
						//  Spiking projections require only post-synaptic component to exist, and for it to receive spikes
						//    - so they require pre-cell to emit spikes
						//  Electrical projections require both pre and post cells to be coupled to twin components of the same synapse type, which requires Vpeer
						//    - so they require both pre and post cells to expose voltage, in the same dimension as syn.component receives
						// 	Continuous projections require only the general rules to be observed
						//    - have fun, be creative
						
						// For example, a STDP synaptic component could be implemented with just a spiking projection, with one post-synaptic component:
						//   then, post-synaptic component should present Voltage and have explicit spike threshold parameter, or should receive spike from its own component as well (allowing event driven simulation)
						
						// NB assume a physical cell has physical dimensionality all over its extent, and presents voltage and receives current
						// 	also assume spike threshold exists (otherwise no triggering is possible), TODO check Vt existence and complain
						// (if it was a cyborg cell it would have been modelled as something coupled with artificial cells in between, probably)
						
						// NB: for SPIKING connections, there can be a different code path for non cell pre...
						
						if(conn.type == Network::Projection::Connection::SPIKING || conn.type == Network::Projection::Connection::ELECTRICAL){
							Int synapse_type = default_synapse_type;
							if(synapse_type < 0){
								
								// find synapse type here
								auto synName = eConn.attribute("synapse").value();
								synapse_type = synaptic_components.get_id(synName);
								if(synapse_type < 0){
									log.error(eConn, "connection synapse type %s not found", synName);
									return false;
								}
								
								// don't forget to send spikes to sole post-synaptic mechanism when voltage exceeds spikeThresh
								
							}
							conn.synapse = synapse_type;
							const auto &syncomp = synaptic_components.get(synapse_type);
							if( conn.type == Network::Projection::Connection::ELECTRICAL ){
								assert( proj.presyn_type == Network::Projection::PresynType::CELL );
								if( !syncomp.HasVpeer(component_types) ){
									log.error(eConn, "connection should use an electrical synapse (using Vpeer)");
									return false;
								}
								// check the interfaces
								
								// symmetric synapse
								if( !CheckSynapticComponentWithCellTypes( log, eConn,
									syncomp, synaptic_components.getName(synapse_type), 
									cell_type_post, cell_type_name_post,
									*cell_type_pre, cell_type_name_pre
								) ) return false;
								if( !CheckSynapticComponentWithCellTypes( log, eConn,
									syncomp, synaptic_components.getName(synapse_type), 
									*cell_type_pre, cell_type_name_pre,
									cell_type_post, cell_type_name_post
								) ) return false;
							}
							if( conn.type == Network::Projection::Connection::SPIKING ){
								if( !syncomp.HasSpikeIn(component_types) ){
									log.error(eConn, "connection should use a synapstic component with a spike in port");
									return false;
								}
								
								if( proj.presyn_type == Network::Projection::PresynType::CELL ){
									// check the interfaces
									if( !CheckSynapticComponentWithCellTypes( log, eConn,
										syncomp, synaptic_components.getName(synapse_type), 
										cell_type_post, cell_type_name_post,
										*cell_type_pre, cell_type_name_pre
									) ) return false;
								}
								else if( proj.presyn_type == Network::Projection::PresynType::EVENT_SERIES ){
									// the spike port will be checked as a segment anyway
								}
								else{ assert(false); return false; }
							}
							
						}
						else if(conn.type == Network::Projection::Connection::CONTINUOUS){
							assert( proj.presyn_type == Network::Projection::PresynType::CELL );
							auto preName = eConn.attribute("preComponent").value();
							auto &synapse_type_pre = conn.continuous.preComponent = synaptic_components.get_id(preName);
							if(conn.continuous.preComponent < 0){
								log.error(eConn, "presynaptic component type %s not found", preName);
								return false;
							}
							auto postName = eConn.attribute("postComponent").value();
							auto &synapse_type_post = conn.continuous.postComponent = synaptic_components.get_id(postName);
							if(conn.continuous.postComponent < 0){
								log.error(eConn, "postsynaptic component type %s not found", postName);
								return false;
							}
							// check the interfaces
							const auto &syncomp_pre  = synaptic_components.get(synapse_type_pre );
							const auto &syncomp_post = synaptic_components.get(synapse_type_post);
							
							// two-part synapse
							if( !CheckSynapticComponentWithCellTypes( log, eConn,
								syncomp_post, postName, 
								cell_type_post, cell_type_name_post,
								*cell_type_pre, cell_type_name_pre
							) ) return false;
							if( !CheckSynapticComponentWithCellTypes( log, eConn,
								syncomp_pre, preName, 
								*cell_type_pre, cell_type_name_pre,
								cell_type_post, cell_type_name_post
							) ) return false;
							
						}
						// now parse some more, mostly common, attributes
						{
						const char *preCell = "preCell";
						const char *preSegment = "preSegment";
						const char *postCell = "postCell";
						const char *postSegment = "postSegment";
						if( uses_old_format ){
							preCell = "preCellId";
							preSegment = "preSegmentId";
							postCell = "postCellId";
							postSegment = "postSegmentId";
						}
						// parse the pre connection differently in each case
						if( proj.presyn_type == Network::Projection::PresynType::EVENT_SERIES ){
						// TODO use different attribute names?
						if(!ParseConnectionPre_EventSeries(log, eConn, *presynaptic_evr, preCell, preSegment, conn.preCell, conn.preSegment, conn.preFractionAlong )) return false; // for pre 
						}
						else if( proj.presyn_type == Network::Projection::PresynType::CELL ){
						if(!ParseConnectionPreOrPost_Cell(log, eConn, *presynaptic_population , pre_morph , preCell , preSegment , "preFractionAlong" , conn.preCell , conn.preSegment , conn.preFractionAlong  )) return false; // for pre 
						// TODO check Vt dependency here for presynaptic cell
						}
						else{ assert(false); return false; }
						// but the post cell is a cell, for now
						if(!ParseConnectionPreOrPost_Cell(log, eConn, postsynaptic_population, post_morph, postCell, postSegment, "postFractionAlong", conn.postCell, conn.postSegment, conn.postFractionAlong )) return false; // for post
						}
						
						if(uses_weight){
							if( !ParseQuantity<Dimensionless>(log, eConn, "weight", conn.weight) ) return false;
						}
						// NB assume continuous projections present no delay, but be ready to handle it anytime LATER
						if(uses_delay){
							if( !ParseQuantity<Time>(log, eConn, "delay", conn.delay) ) return false;
						}
						// could whine if these properties are defined, though unused, LATER
						// TODO also whine if a component with no spike port specifies delay
						
						proj.connections.add(conn, id);
					}
					else{
						//unknown, ignore
					}
				}
				
				// add projection to Network, yay!
				net.projections.add(proj, name);
			}
			else if(strcmp(eNetEl.name(), "explicitInput") == 0){
				const auto &eInp = eNetEl;
				Network::Input input_instance;
				
				//get component
				auto inpName = eInp.attribute("input").value();
				input_instance.component_type = input_sources.get_id(inpName);
				if(input_instance.component_type < 0){
					log.error(eInp, "input type %s not found", inpName);
					return false;
				}
				
				//and positions
				if( !parseCompartmentTarget(
					log, eInp, morphologies, cell_types, net.populations,
					input_instance.population, input_instance.cell_instance, input_instance.segment, input_instance.fractionAlong
				) ){
					return false;
				}
				
				// and check signature compatibility
				Int cell_type_seq = net.populations.get(input_instance.population).component_cell;
				const InputSource &input_source = input_sources.get(input_instance.component_type); // keep specific inpName, just in case
				
				if( !CheckInputComponentWithCellType( log, eInp, input_source, inpName, cell_type_seq ) ) return false; 
				
				net.inputs.push_back(input_instance); // yay!
			}
			else if(strcmp(eNetEl.name(), "inputList") == 0){
				const auto &eInp = eNetEl;
				
				auto list_name = RequiredNmlId(log, eInp, net.input_lists);
				if(!list_name) return false;
				
				Network::InputList list;
				
				//get component
				auto inpName = eInp.attribute("component").value();
				list.component = input_sources.get_id(inpName);
				if(list.component < 0){
					log.error(eInp, "input type %s not found", inpName);
					return false;
				}
				
				auto popName = eInp.attribute("population").value();
				list.population = net.populations.get_id(popName);
				if(list.population < 0){
					log.error(eInp, "input target population %s not found", popName);
					return false;
				}
				
				// and check signature compatibility
				Int cell_type_seq = net.populations.get(list.population).component_cell;
				const InputSource &input_source = input_sources.get(list.component); // keep specific inpName, just in case
				
				if( !CheckInputComponentWithCellType( log, eInp, input_source, inpName, cell_type_seq ) ) return false; 
				
				// TODO verify weights are handled properly
				for(const auto &eInpEl : eInp.children()){
					if( strcmp(eInpEl.name(), "input") == 0 || strcmp(eInpEl.name(), "inputW") == 0 ){
						
						Network::Input input_instance;
						input_instance.component_type = list.component, input_instance.population = list.population;
						Int instance_id;
						
						if(!(StrToL(eInpEl.attribute("id").value(), instance_id) && instance_id >= 0)){
							log.error(eInpEl, "%s requires a non-negative integer id", eInpEl.name());
							return false;
						}
						
						if( !parseCompartmentTarget(log, eInpEl, morphologies, cell_types, net.populations,
							input_instance.population, input_instance.cell_instance, input_instance.segment, input_instance.fractionAlong, true) ) return false;
						
						if( strcmp(eInpEl.name(), "inputW") == 0 ){
							if( !ParseQuantity<Dimensionless>(log, eInpEl, "weight", input_instance.weight) ) return false;
						}
						
						list.input_instances_seq.add((Int) net.inputs.size(), instance_id);
						net.inputs.push_back(input_instance); // yay!
					}
					else{
						// unknown, ignore
					}
				}
				net.input_lists.add(list, list_name);
			}
			// extensions!
			else if(strcmp(eNetEl.name(), "EdenTimeSeriesReader") == 0){
				const auto &eRea = eNetEl;
				
				auto name = RequiredNmlId(log, eRea, net.data_readers);
				if(!name) return false;
				
				Network::TimeSeriesReader reader;
				
				if(!GetSchemeOrRelpathFromUrl(log, eRea, "href", reader.source_url)) return false;
				// NB validate url on backend, until the url format is standardized LATER (or even better, leave it up to the backend?)
				
				const char *sFormat = RequiredAttribute( log, eRea, "format");
				if(!sFormat) return false;
				reader.data_format = sFormat; // NB validate on backend, until the url format is standardized LATER (or even better, leave it up to the backend?)
				
				if(!RequiredPositiveAttribute(log, eRea, "instances", reader.instances)) return false;
				
				const char *sColName = "InputColumn";
				for(const auto &eReaEl : eRea.children()){
					if( strcmp(eReaEl.name(), sColName) == 0 ){
						
						auto name = RequiredNmlId(log, eReaEl, reader.columns);
						if(!name) return false;
						
						Network::TimeSeriesReader::InputColumn column = {Dimension::Unity(), dimensions.GetNative(Dimension::Unity())};
						if(!ParseDimensionAttribute(log, eReaEl, dimensions, column.dimension)) return false;
						const char *sInUnits = "units";
						if(column.dimension != Dimension::Unity() || *eReaEl.attribute(sInUnits).value()){
							const char *sUnits = RequiredAttribute( log, eReaEl, sInUnits);
							if(!sUnits) return false;
							if(!ValidateGetUnits(LogWithElement(log,eReaEl), dimensions, column.dimension, sUnits, column.units)) return false;
						}
						else{
							column.units = dimensions.GetNative(column.dimension); // by default for dimensionless
						}
						reader.columns.add(column, name); // yay!
					}
					else{
						// unknown, ignore
					}
				}
				
				// one last consistency check
				if(reader.columns.contents.empty()){
					log.error(eRea, "%s must have one or more <%s>s", eRea.name(), sColName);
					return false;
				}
				
				net.data_readers.add(reader, name); // yay!
			}
			else if(strcmp(eNetEl.name(), "EdenEventSetReader") == 0){
				const auto &eRea = eNetEl;
				
				auto name = RequiredNmlId(log, eRea, net.event_readers);
				if(!name) return false;
				
				Network::EventSetReader reader;
				
				if(!GetSchemeOrRelpathFromUrl(log, eRea, "href", reader.source_url)) return false;
				// NB validate url on backend, until the url format is standardized LATER (or even better, leave it up to the backend?)
				
				const char *sFormat = RequiredAttribute(log, eRea, "format");
				if(!sFormat) return false;
				reader.data_format = sFormat; // NB validate on backend, until the url format is standardized LATER (or even better, leave it up to the backend?)
				
				if(!RequiredPositiveAttribute(log, eRea, "instances", reader.instances)) return false;
				
				const char *sColName = "Port";
				for(const auto &eReaEl : eRea.children()){
					if( strcmp(eReaEl.name(), sColName) == 0 ){
						
						auto name = RequiredNmlId(log, eReaEl, reader.ports, eRea.name());
						if(!name) return false;
						
						Network::EventSetReader::Port port = {};
						reader.ports.add(port, name); // yay!
					}
					else{
						// unknown, ignore
					}
				}
				
				// one last consistency check
				if(reader.ports.contents.empty()){
					log.error(eRea, "%s must have one or more <%s>s", eRea.name(), sColName);
					return false;
				}
				
				net.event_readers.add(reader, name); // yay!
			}
			else{
				// unknown, ignore
			}
			// printf("done %s\n", eNetEl.name());
		}
		
		networks.add(net, name);
		
		return true;
	}
	
	bool ParseSynapticComponent(const ImportLogger &log, const pugi::xml_node &eSyn){
		
		SynapticComponent syn;
		
		auto name = RequiredNmlId(log, eSyn);
		if(!name) return false;
		
		const char *sType = TagNameOrType(eSyn);
		
		const static NameMap<SynapticComponent::Type> syncomp_types = {
			{"silentSynapse"		  , SynapticComponent::SILENT        },
			{"alphaCurrentSynapse"	  , SynapticComponent::ALPHA_CURRENT },
			{"alphaSynapse"			  , SynapticComponent::ALPHA         },
			{"expOneSynapse"		  , SynapticComponent::EXP           },
			{"expTwoSynapse"		  , SynapticComponent::EXPTWO        },
			{"expThreeSynapse"		  , SynapticComponent::EXPTHREE      },
			
			{"gapJunction"			  , SynapticComponent::GAP           },
			{"linearGradedSynapse"	  , SynapticComponent::GAPHALF       },
			{"gradedSynapse"	      , SynapticComponent::GRADED_PRINZ_ET_AL },
			{"blockingPlasticSynapse" , SynapticComponent::BLOCKING_PLASTIC },
			
			{"expCondSynapse"         , SynapticComponent::PYNN_EXP_COND   },
			{"expCurrSynapse"         , SynapticComponent::PYNN_EXP_CURR },
			{"alphaCondSynapse"       , SynapticComponent::PYNN_ALPHA_COND   },
			{"alphaCurrSynapse"       , SynapticComponent::PYNN_ALPHA_CURR },
		};
		
		// first check if it's one of the channel distribution types
		auto syntype_it = syncomp_types.find( sType );
		if(syntype_it != syncomp_types.end()){
			
			syn.type = syntype_it->second;
			
			if(syn.type == SynapticComponent::SILENT){
				// nothing to do for dummy synapse
			}
			else if(syn.type == SynapticComponent::ALPHA){
				if( !ParseQuantity<Conductance>(log, eSyn, "gbase", syn.alpha.gbase) ) return false;
				if( !ParseQuantity<Voltage>(log, eSyn, "erev", syn.alpha.erev) ) return false;
				if( !ParseQuantity<Time>(log, eSyn, "tau", syn.alpha.tau) ) return false;
			}
			else if(syn.type == SynapticComponent::EXP || syn.type == SynapticComponent::EXPTWO || syn.type == SynapticComponent::EXPTHREE){
				if( !ParseQuantity<Voltage>(log, eSyn, "erev", syn.exp.erev) ) return false;
				
				// tau1, gbase1
				if(syn.type == SynapticComponent::EXP || syn.type == SynapticComponent::EXPTWO){
					if( !ParseQuantity<Time>(log, eSyn, "tauDecay", syn.exp.tauDecay) ) return false;
					if( !ParseQuantity<Conductance>(log, eSyn, "gbase", syn.exp.gbase) ) return false;
				}
				else if(syn.type == SynapticComponent::EXPTHREE){
					if( !ParseQuantity<Time>(log, eSyn, "tauDecay1", syn.exp.tauDecay) ) return false;
					if( !ParseQuantity<Conductance>(log, eSyn, "gbase1", syn.exp.gbase) ) return false;
				}
				
				// tau2
				if(syn.type == SynapticComponent::EXPTWO || syn.type == SynapticComponent::EXPTHREE){
					if( !ParseQuantity<Time>(log, eSyn, "tauRise", syn.exp.tauRise) ) return false;
				}
				
				// tau3
				if( syn.type == SynapticComponent::EXPTHREE){
					if( !ParseQuantity<Time>(log, eSyn, "tauDecay2", syn.exp.tauDecay2) ) return false;
					if( !ParseQuantity<Conductance>(log, eSyn, "gbase2", syn.exp.gbase2) ) return false;
				}
				
				// FIXME complain about tauRise = tauDecay (and support nonetheless)
				
			}
			else if(syn.type == SynapticComponent::GAP || syn.type == SynapticComponent::GAPHALF){
				if( !ParseQuantity<Conductance>(log, eSyn, "conductance", syn.gap.conductance) ) return false;
			}
			else if(syn.type == SynapticComponent::ALPHA_CURRENT){
				if( !ParseQuantity<Current>(log, eSyn, "ibase", syn.alpha_current.ibase) ) return false;
				if( !ParseQuantity<Time>(log, eSyn, "tau", syn.alpha_current.tau) ) return false;
			}
			else if(syn.type == SynapticComponent::GRADED_PRINZ_ET_AL){
				if( !ParseQuantity<Conductance>(log, eSyn, "conductance", syn.graded.gbase ) ) return false;
				if( !ParseQuantity<Voltage>(log, eSyn, "erev"  , syn.graded.erev  ) ) return false;
				if( !ParseQuantity<Voltage>(log, eSyn, "Vth"   , syn.graded.Vth   ) ) return false;
				if( !ParseQuantity<Voltage>(log, eSyn, "delta" , syn.graded.delta ) ) return false;
				if( !ParseQuantity<Frequency>(log, eSyn, "k", syn.graded.k ) ) return false;
			}
			else if(syn.type == SynapticComponent::BLOCKING_PLASTIC){
				if( !ParseQuantity<Conductance>(log, eSyn, "gbase", syn.blopla.gbase) ) return false;
				if( !ParseQuantity<Voltage>(log, eSyn, "erev", syn.blopla.erev) ) return false;
				if( !ParseQuantity<Time>(log, eSyn, "tauDecay", syn.blopla.tauDecay) ) return false;
				if( !ParseQuantity<Time>(log, eSyn, "tauRise", syn.blopla.tauRise) ) return false;
				
				// but also the blocking plastic properties
				// how to tell Blocking LEMS mechanism from Plasticity LEMS mechanism ??
				// perhaps discriminate between base types and hope a child base type is not an ancestor of another base type?
				// it is a mystery
				
				// well, just parse the core types for now
				// TODO check how jLEMS does it
				
				bool has_block_mechanism = false;
				bool has_plasticity_mechanism = false;
				
				// NOTE NeuroMl/NEURON handles event relaying in an ad-hoc way
				//  https://github.com/NeuroML/org.neuroml.export/blob/development/src/main/java/org/neuroml/export/neuron/NeuronWriter.java#L3087
				
				for(auto eSynElm : eSyn.children()){
					
					const char *sElmType = TagNameOrType( eSynElm );
					
					// what is this child element ?
					// perhaps annotation stuff?
					// NOTE wherever there is Standalone, notes, annotation and property tage may exist. Handle properties as they appear, for they must be documented to be used.
					if(
						strcmp( sElmType, "notes" ) == 0
						|| strcmp( sElmType, "annotation" ) == 0
					){
						continue;
					}
					else if( strcmp( sElmType, "voltageConcDepBlockMechanism" ) == 0 ){
						
						if( has_block_mechanism ){
							log.error(eSynElm, "block mechanism already defined" ); // support multiple LATER if the need arises
							return false;
						}
						
						SynapticComponent::BlockingPlasticSynapse::BlockMechanism block;
						block.type = decltype(block)::VOLTAGE_CONC_DEP;
						 
						if( !ParseQuantity<Concentration>(log, eSynElm, "blockConcentration", block.blockConcentration) ) return false;
						if( !ParseQuantity<Concentration>(log, eSynElm, "scalingConc", block.scalingConc) ) return false;
						if( !ParseQuantity<Voltage>(log, eSynElm, "scalingVolt", block.scalingVolt) ) return false;
						
						auto ion_species_name = eSynElm.attribute("species").value();
						if(!*ion_species_name){
							block.species = -1;
							// reject, for now
							log.error(eSynElm, "%s should specify species attribute", sElmType );
							return false;
						}
						else{
							block.species = ion_species.idOrNew(ion_species_name);
						}
						
						syn.blopla.block_mechanism = block;
						has_block_mechanism = true;
					}
					else if(
						strcmp( sElmType, "tsodyksMarkramDepMechanism" ) == 0
						|| strcmp( sElmType, "tsodyksMarkramDepFacMechanism" ) == 0
					){
						if( has_plasticity_mechanism ){
							log.error(eSynElm, "plasticity mechanism already defined" );
							return false;
						}
						
						SynapticComponent::BlockingPlasticSynapse::PlasticityMechanism plasticity;
						
						bool has_fac = strcmp( sElmType, "tsodyksMarkramDepFacMechanism" ) == 0 ;
						if( has_fac ){
							plasticity.type = decltype(plasticity)::TSODYKS_MARKRAM_DEP_FAC;
						}
						else{
							plasticity.type = decltype(plasticity)::TSODYKS_MARKRAM_DEP;
						}
						
						if( !ParseQuantity<Dimensionless>(log, eSynElm, "initReleaseProb", plasticity.initReleaseProb) ) return false;
						if(!( 0 <= plasticity.initReleaseProb && plasticity.initReleaseProb <= 1 )){
							log.error(eSynElm, "initReleaseProb be between 0 and 1" );
							return false;
						}
						if( !ParseQuantity<Time>(log, eSynElm, "tauRec", plasticity.tauRec) ) return false;
						if( has_fac ){
							if( !ParseQuantity<Time>(log, eSynElm, "tauFac", plasticity.tauFac) ) return false;
						}
						
						syn.blopla.plasticity_mechanism = plasticity;
						has_plasticity_mechanism = true;
					}
					else{
						// check for appropriate LEMS components TODO
						log.error(eSynElm, "unknown blocking plastic synapse element %s", sElmType );
						return false;
					}
				}
				
			}
			else if(syn.type == SynapticComponent::PYNN_EXP_COND
				|| syn.type == SynapticComponent::PYNN_ALPHA_COND
			){
				if( !ParseQuantity<Dimensionless>(log, eSyn, "e_rev"  , syn.pynn.e_rev   ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eSyn, "tau_syn", syn.pynn.tau_syn ) ) return false;
			}
			else if(syn.type == SynapticComponent::PYNN_EXP_CURR
				|| syn.type == SynapticComponent::PYNN_ALPHA_CURR
			){
				if( !ParseQuantity<Dimensionless>(log, eSyn, "tau_syn", syn.pynn.tau_syn ) ) return false;
			}
			else{
				log.error(eSyn, "internal error: unknown synapse type %s", sType );
				return false;
			}
			
			// now do a devious trick, refactor LATER
			// along with the LEMS component, attach a mapping to the LEMS version
			// also initialize the id_seq to show it's missing
			syn.component.clear();
			
			if(syn.type == SynapticComponent::SILENT){
				if( !TryLemsifyComponent( log, eSyn, sType, syn.component ) ) return false;
			}
			else if(syn.type == SynapticComponent::GAPHALF){
				if( !TryLemsifyComponent( log, eSyn, sType, {
					{ "conductance", syn.gap.conductance },
				}, syn.component ) ) return false;
			}
			else if(syn.type == SynapticComponent::GAP); // done in Eden core, bother with LATER
			else if(syn.type == SynapticComponent::GRADED_PRINZ_ET_AL){
				if( !TryLemsifyComponent( log, eSyn, sType, {
					{ "erev" , syn.graded.erev  },
					{ "conductance", syn.graded.gbase },
					{ "Vth"  , syn.graded.Vth   },
					{ "delta", syn.graded.delta },
					{ "k"    , syn.graded.k     },
				}, syn.component ) ) return false;
			}
			else if(syn.type == SynapticComponent::ALPHA_CURRENT){
				if( !TryLemsifyComponent( log, eSyn, sType, {
					{ "ibase", syn.alpha_current.ibase },
					{ "tau", syn.alpha_current.tau },
				}, syn.component ) ) return false;
			}
			else if(syn.type == SynapticComponent::ALPHA){
				if( !TryLemsifyComponent( log, eSyn, sType, {
					{ "erev" , syn.alpha.erev  },
					{ "gbase", syn.alpha.gbase },
					{ "tau"  , syn.alpha.tau   },
				}, syn.component ) ) return false;
			}
			else if(syn.type == SynapticComponent::EXP); // done in Eden core, bother with LATER
			else if(syn.type == SynapticComponent::EXPTWO ){
				if( !TryLemsifyComponent( log, eSyn, sType, {
					{ "erev"    , syn.exp.erev     },
					{ "gbase"   , syn.exp.gbase    },
					{ "tauDecay", syn.exp.tauDecay },
					{ "tauRise" , syn.exp.tauRise  },
				}, syn.component ) ) return false;
			}
			else if( syn.type == SynapticComponent::BLOCKING_PLASTIC ){
				if( !TryLemsifyComponent( log, eSyn, "blockingPlasticSynapse", {
					{ "erev"    , syn.blopla.erev     },
					{ "gbase"   , syn.blopla.gbase    },
					{ "tauDecay", syn.blopla.tauDecay },
					{ "tauRise" , syn.blopla.tauRise  },
				}, syn.component ) ) return false;
				
				// also lemsify block mechanism
				if( syn.blopla.block_mechanism.type == SynapticComponent::BlockingPlasticSynapse::BlockMechanism::NONE );
				else if( syn.blopla.block_mechanism.type == SynapticComponent::BlockingPlasticSynapse::BlockMechanism::VOLTAGE_CONC_DEP ){
					auto &block = syn.blopla.block_mechanism;
					if( !TryLemsifyComponent( log, eSyn, "voltageConcDepBlockMechanism", {
						{ "blockConcentration" , block.blockConcentration },
						{ "scalingConc"        , block.scalingConc        },
						{ "scalingVolt"        , block.scalingVolt        },
					}, block.component ) ) return false;
				}
				else{
					log.error(eSyn, "internal error: unknown lemsify block mechanism" );
					return false;
				}
				
				// also lemsify plasticity mechanism
				if( syn.blopla.plasticity_mechanism.type == SynapticComponent::BlockingPlasticSynapse::PlasticityMechanism::NONE );
				else if( syn.blopla.plasticity_mechanism.type == SynapticComponent::BlockingPlasticSynapse::PlasticityMechanism::TSODYKS_MARKRAM_DEP ){
					auto &plasticity = syn.blopla.plasticity_mechanism;
					if( !TryLemsifyComponent( log, eSyn, "tsodyksMarkramDepMechanism", {
						{ "initReleaseProb" , plasticity.initReleaseProb },
						{ "tauRec"          , plasticity.tauRec          },
					}, plasticity.component ) ) return false;
				}
				else if( syn.blopla.plasticity_mechanism.type == SynapticComponent::BlockingPlasticSynapse::PlasticityMechanism::TSODYKS_MARKRAM_DEP_FAC ){
					auto &plasticity = syn.blopla.plasticity_mechanism;
					if( !TryLemsifyComponent( log, eSyn, "tsodyksMarkramDepFacMechanism", {
						{ "initReleaseProb" , plasticity.initReleaseProb },
						{ "tauRec"          , plasticity.tauRec          },
						{ "tauFac"          , plasticity.tauFac          },
					}, plasticity.component ) ) return false;
				}
				else{
					log.error(eSyn, "internal error: unknown lemsify plasticity mechanism" );
					return false;
				}
				
			}
			else if(syn.type == SynapticComponent::EXPTHREE){
				if( !TryLemsifyComponent( log, eSyn, sType, {
					{ "erev"     , syn.exp.erev      },
					{ "gbase1"   , syn.exp.gbase     },
					{ "gbase2"   , syn.exp.gbase2    },
					{ "tauDecay1", syn.exp.tauDecay  },
					{ "tauRise"  , syn.exp.tauRise   },
					{ "tauDecay2", syn.exp.tauDecay2 },
				}, syn.component ) ) return false;
			}
			else if(
				syn.type == SynapticComponent::PYNN_EXP_COND
				|| syn.type == SynapticComponent::PYNN_ALPHA_COND
			){
				if( !TryLemsifyComponent( log, eSyn, sType, {
					{ "e_rev"  , syn.pynn.e_rev   },
					{ "tau_syn", syn.pynn.tau_syn },
				}, syn.component ) ) return false;
			}
			else if(
				syn.type == SynapticComponent::PYNN_EXP_CURR
				|| syn.type == SynapticComponent::PYNN_ALPHA_CURR
			){
				if( !TryLemsifyComponent( log, eSyn, sType, {
					{ "tau_syn", syn.pynn.tau_syn },
				}, syn.component ) ) return false;
			}
			else{
				log.error(eSyn, "internal error: unknown lemsify synaptic component" );
				return false;
			}
			
		}
		else{
			// should be a LEMS component
			// perhaps check parm completion for the LEMSified core components too? nah, let the NANs leak instead (it would be an internal bug anyway)
			
			if( !ParseComponentInstanceSynapticComponent(log, eSyn, component_types, dimensions, sType, syn.component) ) return false;
			
				// printf("Tagifieddd %s\n",sType);
			syn.type = decltype(syn)::COMPONENT;
			// log.error(eSyn, "unknown synapse type %s", eSyn.name() );
			// return false;
		}
		
		//add synapse to synaptic components! yay!
		synaptic_components.add(syn, name);
		
		return true;
	}
	
	
	bool ValidateComponentInstanceIndependentArtificialCell( const ImportLogger &log, const pugi::xml_node &eReferred, const ComponentInstance &instance ) const {
		
		std::map< std::string, ComponentType::Requirement > provided_requirements;
		CoverCommonGlobalRequirements( provided_requirements );
		CoverCommonRequirement( "iSyn", LEMS_Current, provided_requirements );
		CoverCommonRequirement( "ISyn", Dimension::Unity(), provided_requirements );
		
		std::map< std::string, ComponentType::Requirement > required_exposures   ;
		// situational requirements will be resolved at synaptic/input projection time
		
		std::map< std::string, ComponentType::EventPortIn > provided_event_inputs;
		CoverCommonEventIn( "in", provided_event_inputs ); // just in case, no harm in adding a placeholder
		
		std::map< std::string, ComponentType::EventPortOut > required_event_outputs;
		// situational requirements will be resolved at synaptic/input projection time
		
		const auto &comp = component_types.get(instance.id_seq);
		const char *type = component_types.getName(instance.id_seq);
		
		if( !ValidateComponentTypeInterface(log, eReferred, comp, dimensions, type, provided_requirements, required_exposures, provided_event_inputs, required_event_outputs ) ) return false;
		if( !ValidateComponentInstanceCompleteness(log, eReferred, comp, type, instance) ) return false;
		
		return true;
	}
	
	bool RealizeComponentAsArtificialCell(const ImportLogger &log, const pugi::xml_node &eRealizationPoint, const ComponentInstance &compinst, const char *name ){
		if( !ValidateComponentInstanceIndependentArtificialCell( log, eRealizationPoint, compinst ) ) return false;
		
		// perhaps check interface more, LATER
		
		// it can be an artificial cell, proceed
		if( !cell_types.has(name) ){
			CellType new_cell_type;
			new_cell_type.type = CellType::ARTIFICIAL;
			new_cell_type.artificial.type = ArtificialCell::COMPONENT;
			new_cell_type.artificial.component = compinst;
			
			cell_types.add( new_cell_type, name );
		}
		else{
			// perhaps check if the same?
			// though this shouldn't happen if this routine was called, due to missing a cell type with that name
			log.error( eRealizationPoint, "is already an existing cell type" );
			return false;
			
		}
		
		return true; // yay !
	}
	
	bool ValidateComponentInstanceIndependentSpikeSource( const ImportLogger &log, const pugi::xml_node &eReferred, const ComponentInstance &instance ) const {
		
		std::map< std::string, ComponentType::Requirement > provided_requirements;
		CoverCommonGlobalRequirements( provided_requirements );
		
		std::map< std::string, ComponentType::Requirement > required_exposures   ;
		//CoverCommonRequirement( "i", LEMS_Current, required_exposures );
		
		std::map< std::string, ComponentType::EventPortIn > provided_event_inputs;
		
		std::map< std::string, ComponentType::EventPortOut > required_event_outputs;
		CoverCommonEventOut( "spike", required_event_outputs );
		
		const auto &comp = component_types.get(instance.id_seq);
		const char *type = component_types.getName(instance.id_seq);
		
		if( !ValidateComponentTypeInterface(log, eReferred, comp, dimensions, type, provided_requirements, required_exposures, provided_event_inputs, required_event_outputs ) ) return false;
		if( !ValidateComponentInstanceCompleteness(log, eReferred, comp, type, instance) ) return false;
		
		return true;
	}
	// LATER add a 'Try' version of this
	bool RealizeInputSourceAsArtificialCell(const ImportLogger &log, const pugi::xml_node &eRealizationPoint, const InputSource &inp, const char *name ){
		
		Dimension unused;
		if( inp.type == InputSource::COMPONENT ){
			
			if( !ValidateComponentInstanceIndependentSpikeSource( log, eRealizationPoint, inp.component ) ) return false;
			
			// fall through to cell type creation
		}
		else{
			// core type
			
			if( inp.GetCurrentOutputAndDimension( component_types, synaptic_components, unused ) ){
				log.error( eRealizationPoint, "input type %s generates current, so it cannot be an individual artificial cell", name );
				return false;
			}
			else if(
				inp.HasSpikeOut( component_types )
			){
				// fall through to cell type creation
			} 
			else{
				// not covered, catch-all case
				log.error( eRealizationPoint, "input type %s cannot be an individual artificial cell", name );
				return false;
			}
		}
		
		// it can be an artificial cell, proceed
		if( !cell_types.has(name) ){
			CellType new_cell_type;
			new_cell_type.type = CellType::ARTIFICIAL;
			new_cell_type.artificial.type = ArtificialCell::SPIKE_SOURCE;
			new_cell_type.artificial.spike_source_seq = input_sources.get_id(name);
			
			cell_types.add( new_cell_type, name );
		}
		else{
			// perhaps check if the same?
			// though this shouldn't happen if this routine was called, due to missing a cell type with that name
			log.error( eRealizationPoint, "is already an existing cell type" );
			return false;
			
		}
		
		return true; // yay !
	}
	bool ParseInputSource(const ImportLogger &log, const pugi::xml_node &eInp){
		
		InputSource inp;
		
		auto name = RequiredNmlId(log, eInp);
		if(!name) return false;
		
		auto sType = TagNameOrType(eInp);
		
		const static NameMap<InputSource::Type> input_types = {
			{"pulseGenerator"     , InputSource::PULSE                },
			{"pulseGeneratorDL"   , InputSource::PULSE_DL             },
			{"sineGenerator"      , InputSource::SINE                 },
			{"sineGeneratorDL"    , InputSource::SINE_DL              },
			{"rampGenerator"      , InputSource::RAMP                 },
			{"rampGeneratorDL"    , InputSource::RAMP_DL              },
			{"voltageClamp"       , InputSource::VOLTAGE_CLAMP        },
			{"voltageClampTriple" , InputSource::VOLTAGE_CLAMP_TRIPLE },
			{"timedSynapticInput" , InputSource::TIMED_SYNAPTIC       },
			{"poissonFiringSynapse" , InputSource::POISSON_SYNAPSE      },
			{"transientPoissonFiringSynapse" , InputSource::POISSON_SYNAPSE_TRANSIENT },

			{"spikeArray"         , InputSource::SPIKE_LIST           },
			{"spikeGenerator"     , InputSource::SPIKE_PERIODIC       },
			{"spikeGeneratorRandom" , InputSource::SPIKE_RANDOM       },
			{"spikeGeneratorPoisson" , InputSource::SPIKE_POISSON       },
			{"spikeGeneratorRefPoisson" , InputSource::SPIKE_POISSON_REF  },
			{"SpikeSourcePoisson" , InputSource::PYNN_SPIKE_POISSON       },
		};
		
		auto GetSpikesForElement = [  ]( const ImportLogger &log, const pugi::xml_node &eInp, std::vector<InputSource::Spike> &spikes ){
			
			// gather a list of spikes, at the time they occur
			for(auto eInpEl : eInp.children()){
				if( strcmp(eInpEl.name(), "spike") == 0 ){
					InputSource::Spike spike;
					if( !ParseQuantity<Time>(log, eInpEl, "time", spike.time_of_occurrence ) ) return false;
					spikes.push_back(spike);
				}
				else{
					// unknown, ignore
					// incl. Standalone sub-elements
				}
			}
			
			// sort the spikes by time after they are input
			std::sort(spikes.begin(), spikes.end(), [](auto a, auto b){
				return a.time_of_occurrence < b.time_of_occurrence;
			});
			
			return true;
		};
		auto GetSynapticSubcomponent = [ this ]( const ImportLogger &log, const pugi::xml_node &eInp, Int &syn_seq){
			
			const char *sSynapse = RequiredAttribute( log, eInp, "synapse");
			if(!sSynapse) return false;
			syn_seq = synaptic_components.get_id(sSynapse);
			if( syn_seq < 0 ){
				log.error(eInp, "unknown synapse type %s", sSynapse );
				return false;
			}
			// verify spikeTarget path LATER, who cares anyway
			const auto &syn = synaptic_components.get(syn_seq);
			
			// also check interface
			if( !syn.HasSpikeIn(component_types) ){
				log.error(eInp, "synapse should receive spike events", sSynapse );
				return false;
			}
			if( syn.HasVpeer(component_types) ){
				log.error(eInp, "synapse should not depend on peer voltage", sSynapse );
				return false;
			}
			
			return true;
		};
		
		// first check if it's one of the channel distribution types
		auto inputype_it = input_types.find( sType );
		if(inputype_it != input_types.end()){
			inp.type = inputype_it->second;
			
			if(inp.type == InputSource::PULSE){
				if( !ParseQuantity<Current      >(log, eInp, "amplitude", inp.amplitude ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "delay"    , inp.delay     ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "duration" , inp.duration  ) ) return false;
			}
			else if(inp.type == InputSource::PULSE_DL){
				if( !ParseQuantity<Dimensionless>(log, eInp, "amplitude", inp.amplitude ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "delay"    , inp.delay     ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "duration" , inp.duration  ) ) return false;
			}
			else if(inp.type == InputSource::SINE){
				if( !ParseQuantity<Current      >(log, eInp, "amplitude", inp.amplitude ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "delay"    , inp.delay     ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "duration" , inp.duration  ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "period"   , inp.period    ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eInp, "phase"    , inp.phase     ) ) return false;
			}
			else if(inp.type == InputSource::SINE_DL){
				if( !ParseQuantity<Dimensionless>(log, eInp, "amplitude", inp.amplitude ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "delay"    , inp.delay     ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "duration" , inp.duration  ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "period"   , inp.period    ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eInp, "phase"    , inp.phase     ) ) return false;
			}
			else if(inp.type == InputSource::RAMP){
				if( !ParseQuantity<Current      >(log, eInp, "baselineAmplitude", inp.amplitude       ) ) return false;
				if( !ParseQuantity<Current      >(log, eInp, "startAmplitude"   , inp.startAmplitude  ) ) return false;
				if( !ParseQuantity<Current      >(log, eInp, "finishAmplitude"  , inp.finishAmplitude ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "delay"            , inp.delay           ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "duration"         , inp.duration        ) ) return false;
			}
			else if(inp.type == InputSource::RAMP_DL){
				if( !ParseQuantity<Dimensionless>(log, eInp, "baselineAmplitude", inp.amplitude       ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eInp, "startAmplitude"   , inp.startAmplitude  ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eInp, "finishAmplitude"  , inp.finishAmplitude ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "delay"            , inp.delay           ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "duration"         , inp.duration        ) ) return false;
			}
			else if(inp.type == InputSource::VOLTAGE_CLAMP){
				if( !ParseQuantity<Voltage      >(log, eInp, "targetVoltage", inp.testingVoltage ) ) return false;
				
				// an unfortunate incident means seriesResistance is off-scale on simple voltageClamp
				if ( eInp.attribute("simpleSeriesResistance") ){
					if( !ParseQuantity<Resistance   >(log, eInp, "simpleSeriesResistance"   , inp.seriesResistance ) ) return false;
				}
				else if( eInp.attribute("seriesResistance") ){
					if( !ParseQuantity<Resistance   >(log, eInp, "seriesResistance"   , inp.seriesResistance ) ) return false;
					inp.seriesResistance *= 100000; // it happens
				}
				else{
					log.error(eInp, "attribute simpleSeriesResistance or seriesResistance required");
					return false;
				}
				
				if( !ParseQuantity<Time         >(log, eInp, "delay"    , inp.delay     ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "duration" , inp.duration  ) ) return false;
			}
			else if(inp.type == InputSource::VOLTAGE_CLAMP_TRIPLE){
				if( !ParseQuantity<Time         >(log, eInp, "delay"                  , inp.delay               ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "duration"               , inp.duration            ) ) return false;
				if( !ParseQuantity<Dimensionless>(log, eInp, "active"                 , inp.active              ) ) return false;
				if( !ParseQuantity<Voltage      >(log, eInp, "conditioningVoltage"    , inp.conditioningVoltage ) ) return false;
				if( !ParseQuantity<Voltage      >(log, eInp, "testingVoltage"         , inp.testingVoltage      ) ) return false;
				if( !ParseQuantity<Voltage      >(log, eInp, "returnVoltage"          , inp.returnVoltage       ) ) return false;
				if( !ParseQuantity<Resistance   >(log, eInp, "simpleSeriesResistance" , inp.seriesResistance    ) ) return false;
			}
			else if(inp.type == InputSource::TIMED_SYNAPTIC){
				
				// get spikes
				if( !GetSpikesForElement( log, eInp, inp.spikes ) ) return false;
				
				// and attached synapse
				if( !GetSynapticSubcomponent( log, eInp, inp.synapse ) ) return false;
				
				// that's all
			}
			else if(inp.type == InputSource::POISSON_SYNAPSE){
				if( !ParseQuantity<Frequency    >(log, eInp, "averageRate", inp.averageRate ) ) return false;
				
				// and attached synapse
				if( !GetSynapticSubcomponent( log, eInp, inp.synapse ) ) return false;
			}
			else if(inp.type == InputSource::POISSON_SYNAPSE_TRANSIENT){
				if( !ParseQuantity<Frequency    >(log, eInp, "averageRate", inp.averageRate ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "delay"      , inp.delay       ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "duration"   , inp.duration    ) ) return false;
				
				// and attached synapse
				if( !GetSynapticSubcomponent( log, eInp, inp.synapse ) ) return false;
			}
			else if(inp.type == InputSource::SPIKE_LIST){
				if( !GetSpikesForElement( log, eInp, inp.spikes ) ) return false;
			}
			else if(inp.type == InputSource::SPIKE_PERIODIC){
				if( !ParseQuantity<Time         >(log, eInp, "period"   , inp.period    ) ) return false;
			}
			else if(inp.type == InputSource::SPIKE_RANDOM){
				if( !ParseQuantity<Time         >(log, eInp, "maxISI"   , inp.maxISI    ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "minISI"   , inp.minISI    ) ) return false;
			}
			else if(inp.type == InputSource::SPIKE_POISSON){
				if( !ParseQuantity<Frequency    >(log, eInp, "averageRate", inp.averageRate    ) ) return false;
			}
			else if(inp.type == InputSource::SPIKE_POISSON_REF){
				if( !ParseQuantity<Frequency    >(log, eInp, "averageRate", inp.averageRate ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "minimumISI" , inp.minISI      ) ) return false;				
			}
			else if(inp.type == InputSource::PYNN_SPIKE_POISSON){
				if( !ParseQuantity<Frequency    >(log, eInp, "rate"     , inp.averageRate    ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "start"    , inp.delay     ) ) return false;
				if( !ParseQuantity<Time         >(log, eInp, "duration" , inp.duration  ) ) return false;
			}
			else{
				log.error(eInp, "internal error: unknown input type %s", sType );
				return false;
			}
			
			// now do a devious trick, refactor much LATER if ever
			// along with the LEMS component, attach a mapping to the LEMS version
			// also initialize the id_seq to show it's missing
			inp.component.clear();
			if(inp.type == InputSource::PULSE); // done in Eden core, bother with LATER
			else if(inp.type == InputSource::PULSE_DL){
				if( !TryLemsifyComponent( log, eInp, sType, {
					{ "amplitude", inp.amplitude },
					{ "duration" , inp.duration  },
					{ "delay"    , inp.delay     },
				}, inp.component ) ) return false;
			}
			else if(inp.type == InputSource::SINE
				|| inp.type == InputSource::SINE_DL
			){
				if( !TryLemsifyComponent( log, eInp, sType, {
					{ "amplitude", inp.amplitude },
					{ "duration" , inp.duration  },
					{ "delay"    , inp.delay     },
					{ "period"   , inp.period    },
					{ "phase"    , inp.phase     },
				}, inp.component ) ) return false;
			}
			else if(inp.type == InputSource::RAMP
				|| inp.type == InputSource::RAMP_DL
			){
				if( !TryLemsifyComponent( log, eInp, sType, {
					{ "baselineAmplitude", inp.amplitude       },
					{ "duration"         , inp.duration        },
					{ "delay"            , inp.delay           },
					{ "startAmplitude"   , inp.startAmplitude  },
					{ "finishAmplitude"  , inp.finishAmplitude },
				}, inp.component ) ) return false;
			}
			else if(inp.type == InputSource::VOLTAGE_CLAMP){
				if( !TryLemsifyComponent( log, eInp, sType, {
					{ "duration"               , inp.duration         },
					{ "delay"                  , inp.delay            },
					{ "targetVoltage"          , inp.testingVoltage   },
					{ "simpleSeriesResistance" , inp.seriesResistance },
				}, inp.component ) ) return false;
			}
			else if(inp.type == InputSource::VOLTAGE_CLAMP_TRIPLE){
				if( !TryLemsifyComponent( log, eInp, sType, {
					{ "duration"               , inp.duration            },
					{ "delay"                  , inp.delay               },
					{ "active"                 , inp.active              },
					{ "conditioningVoltage"    , inp.conditioningVoltage },
					{ "testingVoltage"         , inp.testingVoltage      },
					{ "returnVoltage"          , inp.returnVoltage       },
					{ "simpleSeriesResistance" , inp.seriesResistance    },
				}, inp.component ) ) return false;
			}
			else if(inp.type == InputSource::POISSON_SYNAPSE){
				if( !TryLemsifyComponent( log, eInp, sType, {
					{ "averageRate", inp.averageRate },
				}, inp.component ) ) return false;
			}
			else if(inp.type == InputSource::POISSON_SYNAPSE_TRANSIENT){
				if( !TryLemsifyComponent( log, eInp, sType, {
					{ "averageRate", inp.averageRate },
					{ "duration"   , inp.duration    },
					{ "delay"      , inp.delay       },
				}, inp.component ) ) return false;
			}
			else if(inp.type == InputSource::SPIKE_PERIODIC){
				if( !TryLemsifyComponent( log, eInp, sType, {
					{ "period"               , inp.period         },
				}, inp.component ) ) return false;
			}
			else if(inp.type == InputSource::SPIKE_RANDOM){
				if( !TryLemsifyComponent( log, eInp, sType, {
					{ "maxISI"               , inp.maxISI         },
					{ "minISI"               , inp.minISI         },
				}, inp.component ) ) return false;
			}
			else if(inp.type == InputSource::SPIKE_POISSON){
				if( !TryLemsifyComponent( log, eInp, sType, {
					{ "averageRate"          , inp.averageRate    },
				}, inp.component ) ) return false;
			}
			else if(inp.type == InputSource::SPIKE_POISSON_REF){
				if( !TryLemsifyComponent( log, eInp, sType, {
					{ "averageRate"          , inp.averageRate    },
					{ "minimumISI"               , inp.minISI         },
				}, inp.component ) ) return false;
			}
			else if(inp.type == InputSource::PYNN_SPIKE_POISSON){
				if( !TryLemsifyComponent( log, eInp, sType, {
					{ "rate"          , inp.averageRate },
					{ "start"         , inp.delay       },
					{ "duration"      , inp.duration    },
				}, inp.component ) ) return false;
			}
			if(inp.type == InputSource::PULSE); // done in Eden core only
			else{
				// pass through, i guess? should be implemented as core component, though
			}
			
		}
		else{
			// perhaps a LEMS component?
			inp.type = InputSource::COMPONENT;
			if( !ParseComponentInstanceInputCurrent(log, eInp, component_types, dimensions, sType, inp.component) ) return false;
		}
		
		//add input type to input sources! yay!
		input_sources.add(inp, name);
		
		// NB a pseudo-input that is a pure spike source (without current, or voltage dependency etc.) can also be considered an artificial cell !
		
		return true;
	}
	bool ParseSimulation_CustomSetupFile(const ImportLogger &log, const pugi::xml_node &eSimEl, const Network &net, const char *sRelFilename, Simulation::CustomSetup &custom_init ){
		// to avoid more complicated coding, load the whole file in memory until it's time for more efficiency (in which case a binary format can be used instead).
		std::vector< std::vector<std::string> > tokenized_lines;
		// TODO tokenized_offsets
		std::vector< int > tokenized_line_number;
		std::vector< ptrdiff_t > tokenized_line_byte_offset;
		
		const char *loading_from_file = log.GetFilenameFromElement(eSimEl);
		std::string filename_s = GetRelativeFilePath((loading_from_file ? loading_from_file : "."), sRelFilename);
		const char *filename = filename_s.c_str();
		
		std::ifstream fin(filename, std::ios::binary); // binary mode, to count CRLF as two bytes
		// check if opening a file failed
		if (fin.fail()) {
			log.error(eSimEl, "could not open file \"%s\": %s", filename, strerror(errno));
			fin.close();
			return false;
		}
		int source_line = 1;
		ptrdiff_t byte_offset = 0;
		for( std::string line; std::getline(fin, line, '\n'); source_line++, byte_offset += line.size()+1){
			// printf("line: %s\n", line.c_str());
			// trim whitespace, to be nice; also multiple whitespaces in a row, to reduce friction
			// nice eh? https://stackoverflow.com/questions/2275160/splitting-a-string-by-whitespace-in-c#comment118452224_2275160
			std::istringstream streamline(line);
			std::vector<std::string> tokens{std::istream_iterator<std::string>(streamline), {}};
			// for (auto token : tokens ){ printf(" %s", token.c_str()); } puts("\n");
			if( tokens.empty() ) continue; // empty line
			if( !tokens[0].empty() && tokens[0][0] == '#') continue; // comment line
			tokenized_lines.push_back(tokens);
			tokenized_line_number.push_back(source_line);
			tokenized_line_byte_offset.push_back(byte_offset);
		}
		// close the file stream
		fin.close();
		
		auto &statements = custom_init.statements;
		
		{
		const auto &that_log = log;
		
		struct LogCustomFile{
			const ImportLogger &log;
			const char *filename;
			const std::vector< ptrdiff_t > &tokenized_line_byte_offset;
			LogCustomFile(const ImportLogger &_l, const char *_f, const std::vector< ptrdiff_t > &_lo): log(_l), filename(_f), tokenized_line_byte_offset(_lo) {}
			
			void warning(Int tokenized_lineno, Int token, const char *format, ...) const {
				va_list args; va_start(args, format);
				log.warning(args, filename, tokenized_line_byte_offset[tokenized_lineno], format);
				va_end(args);
			}
			void error(Int tokenized_lineno, Int token, const char *format, ...) const {
				va_list args; va_start(args, format);
				log.error(args, filename, tokenized_line_byte_offset[tokenized_lineno], format);
				va_end(args);
			}
			void warning(va_list args, Int tokenized_lineno, Int token, const char *format) const {
				log.warning(args, filename, tokenized_line_byte_offset[tokenized_lineno], format);
			}
			void error(va_list args, Int tokenized_lineno, Int token, const char *format) const {
				log.error(args, filename, tokenized_line_byte_offset[tokenized_lineno], format);
			}
		} log = {that_log, filename, tokenized_line_byte_offset};
		
		struct LogInsideToken : public ILogProxy{
			const LogCustomFile &log;
			Int tokenized_lineno, token;
			LogInsideToken(const LogCustomFile &_l, Int _n, Int _t): log(_l), tokenized_lineno(_n), token(_t) {}
			void error(const char *format, ...) const {
				va_list args; va_start(args, format);
				log.error(args, tokenized_lineno, token, format);
				va_end (args);
			}
			void warning(const char *format, ...) const {
				va_list args; va_start(args, format);
				log.warning(args, tokenized_lineno, token, format);
				va_end (args);
			}
		};
		
		// when container is a CollectionWithIds or implements the same interface for getSequential, getId
		auto IdListToSeqList = [&log](Int lineno, Int token, const std::string &targets_list, const auto &container, const char *item_type, const char *container_type, auto &items_seq ){
			auto ttokens = string_split(targets_list, ",");
			for( const std::string &ttoken : ttokens){
				long cell_id = -1, item_seq = -1;
				if(!( StrToL(ttoken.c_str(), cell_id) && (item_seq = container.getSequential(cell_id)) >= 0 )){
					log.error(lineno, token, "could not resolve %s with id \"%s\" in %s", item_type, ttoken.c_str(), container_type);
					return false;
				}
				items_seq.push_back(item_seq);
			}
			// check if the list has duplicates, that would be ambiguous !
			// TODO extract items list
			auto items_sorted = items_seq;
			std::sort(items_sorted.begin(), items_sorted.end());
			for(Int i = 0; i+1 < (Int) items_sorted.size(); i++){
				if(items_sorted[i] == items_sorted[i+1]){
					log.error(lineno, token, "%s id \"%s\" is specified more than once", item_type, std::to_string(container.getId(items_sorted[i])).c_str());
					return false;
				}
			}
			// TODO check if set is empty, i guess?
			return true;
		};
		// when container is a std::vector
		auto SeqListToSeqList = [&IdListToSeqList](Int lineno, Int token, const std::string &targets_list, const auto &container, const char *item_type, const char *container_type, auto &items_seq ){
			struct LikeACollectionWithIds{
				decltype(container) &_container;
				Int getSequential(Int id) const {
					if( 0 <= id && id <= (Int)_container.size()) return id;
					else return -1;
				}
				Int getId(Int id) const { return id; }
			} confaker = {container};
			return IdListToSeqList(lineno, token, targets_list, confaker, item_type, container_type, items_seq);
		};
		(void) SeqListToSeqList; // in case direct-mapped vectors appear LATER
		
		auto ParseCheckRef = [this, &net](const auto &log, const std::string &value_term, const Dimension &dimension, Simulation::LemsQuantityPath &path_value ){
						
			if(!model.ParseLemsQuantityPath(log, value_term.c_str(), net, path_value)) return false;
			// check type and dimension
			ComponentType::NamespaceThing::Type path_value_type; Dimension value_dimension;
			if(!model.GetLemsQuantityPathType(net, path_value, path_value_type, value_dimension)) return false;
			if(path_value_type != ComponentType::NamespaceThing::STATE){
				log.error("value for VariableReference must refer to a state variable");
			}
			if( value_dimension != dimension ){
				log.error("value for VariableReference has dimension %s, but it should have dimension %s like the property being set", dimensions.Stringify(value_dimension), dimensions.Stringify(dimension));
			}
			return true;
		};
		
		auto MustBeSettable = [](const auto &log, ComponentType::NamespaceThing::Type &path_type, bool &is_ref){
			is_ref = false;
			if(path_type == ComponentType::NamespaceThing::STATE
			|| path_type == ComponentType::NamespaceThing::PROPERTY){
				// it's ok, they can be set
			}
			else if(path_type == ComponentType::NamespaceThing::VARREQ){
				is_ref = true;
			}
			else{
				log.error("term must refer to a State variable, Parameter, Property, or VariableRequirement; not a %s", ComponentType::NamespaceThing::getTypeName(path_type));
				return false;
			}
			return true;
		};
					
		auto ParseToken_ValueSpecifier = [this, &ParseCheckRef](const auto &log, Int lineno, const std::vector<std::string> &line, const Dimension &dimension, bool is_ref, Simulation::CustomSetup::Statement &set, Real &value, LemsUnit &real_value_units, Simulation::LemsQuantityPath &path_value, Int &tokens_consumed){
			
			Int value_token = tokens_consumed;
			const std::string &value_term = line[value_token];
			tokens_consumed++;
			
			if(value_term == "multi"){
				set.multi_mode = true;
			}
			else if(value_term == "cable"){
				set.cable_mode = true;
			}
			else if(value_term == "multicable"){
				set.multi_mode = true;
				set.cable_mode = true;
			}
			else{
				// it should be a value
				if(is_ref){
					LogInsideToken log_proxy = {log, lineno, value_token};
					if(!ParseCheckRef(log_proxy, value_term, dimension, path_value)) return false;
				}
				else{
					if(!StrToF(value_term.c_str(), value)){
						log.error(lineno, value_token, "could not parse %s as numerical value, 'multi', 'cable' or 'cablemulti'", value_term.c_str());
					}
				}
			}
			
			// then if the value is dimensional, get units
			if(	!is_ref && dimension != Dimension::Unity() ){
				Int units_token = tokens_consumed;
				if(!( units_token < (Int)line.size() )){
					log.error(lineno, value_token, "statement needs units for property of dimensionality %s", dimensions.Stringify(dimension).c_str());
					return false;
				}
				const std::string &unit_name = line[units_token];
				tokens_consumed++;
				
				// see also ParseLemsQuantity, LATER refer to this logic and perhaps extract it for when validating supplied units
				if(!dimensions.Has(dimension)){
					log.error(lineno, units_token, "there are no specified %s units", dimensions.Stringify(dimension).c_str() );
					return false;
				}
				bool units_ok = false; // TODO ValidateGetUnits
				for( auto scale : dimensions.GetUnits(dimension) ){
					if(unit_name == scale.name){
						real_value_units = scale;
						value = scale.ConvertTo( value, dimensions.GetNative(dimension) ); // adjust units right away, if it's a single value.  If it's multi, adjust as the values rows are being read. 
						units_ok = true;
						break;
					}
				}
				if(!units_ok){
					//unit name not found in list!
					std::string known_list; for( auto scale : dimensions.GetUnits(dimension) ) known_list += " ", known_list += scale.name;
					log.error(lineno, units_token, "unknown property units: %s for %s (supported:%s)", unit_name.c_str(), dimensions.Stringify(dimension).c_str(), known_list.c_str() );
					return false;
				}
			}
			
			// done
			return true;
		};
		
		// read this many lines (or use value or path_value if inline) into the set statement, according to its data flags
		auto ParseSetData = [&ParseCheckRef, &tokenized_lines, this]( 
			const auto &log, Int lineno, const Int all_instances_count, 
			const Real &value, const Simulation::LemsQuantityPath &path_value, const ScaleEntry &real_value_units, const Dimension &dimension, const bool &is_ref, 
			Simulation::CustomSetup::Statement &set, Int &lines_to_advance 
		){
			// then if multi was specified, read this many lines of data
			auto ReadDataRow_Number = [&log]( Int linenono, Int columns_to_read, const auto &line, std::vector<Real> &row){
				assert((Int)line.size() == 1+columns_to_read);
				
				row.resize(columns_to_read);
				for(Int i = 0; i < columns_to_read; i++){
					const auto &ttoken = line[1+i];
					if(!StrToF(ttoken.c_str(), row[i])){
						log.error(linenono, 1+i, " could not parse as number: %s", ttoken.c_str() );
						return false;
					}
				}
				return true;
			};
			auto ReadDataRow_Ref = [&log, &ParseCheckRef]( Int linenono, Int columns_to_read, const Dimension &dimension, const auto &line, std::vector<Simulation::LemsQuantityPath> &row){
				assert((Int)line.size() == 1+columns_to_read);
				
				row.resize(columns_to_read);
				for(Int i = 0; i < columns_to_read; i++){
					const auto &ttoken = line[1+i];
					LogInsideToken log_proxy = {log, linenono, 1+i};
					if(!ParseCheckRef(log_proxy, ttoken, dimension, row[i])) return false;
				}
				return true;
			};
			
			Int items_to_read = 0;
			if(set.multi_mode){
				if(set.items_seq.empty()) items_to_read = all_instances_count;
				else items_to_read = (Int) set.items_seq.size();
			}
			
			Int lines_to_read = 0, columns_to_read = -1;
			if(set.cable_mode){
				columns_to_read = -1;
				
				if(set.multi_mode) lines_to_read = 1+items_to_read;
				else lines_to_read = 1+1; // one for fractionAlong cable, one for values
			}
			else if(set.multi_mode){
				lines_to_read = 1;
				columns_to_read = items_to_read;
			}
			else{
				// it must be a single value
				if(is_ref){
					set.ref_data.resize(1);
					set.ref_data[0].resize(1);
					set.ref_data[0][0] = path_value;
				}
				else{
					if( value != value ){
						log.error(lineno, 0, "internal error: single value %s unset", ""); // LATER bother passing the token or sth, because its location may vary
						return false;
					}
					set.real_data.resize(1);
					set.real_data[0].resize(1);
					set.real_data[0][0] = value;
				}
			}
			// printf("columns %ld\n", columns_to_read);
			for(Int linenono = lineno + 1; linenono < lineno + 1 + lines_to_read; linenono++){
				if(!( linenono < (Int) tokenized_lines.size() && !tokenized_lines[linenono].empty() && tokenized_lines[linenono][0] == "values" )){
					log.error(std::min(linenono-1, (Int) (tokenized_lines.size() - 1)), 0, "'set' line has only %lld 'values' lines following, when it should have %lld", (long long) (linenono - (lineno + 1)), (long long)lines_to_read);
					return false;
				}
				
				const auto &tokens = tokenized_lines[linenono];
				assert(tokens.size() > 0);
				
				std::vector<Real> new_row;
				
				// Whatever the type of value, in cable mode the first row must be real data of fractionAlong
				if(set.cable_mode && linenono == lineno){
					// the following rows will have as many elements as this index row.
					columns_to_read = (Int)tokens.size()-1;
					
					if(!ReadDataRow_Number(linenono, columns_to_read, tokens, new_row )) return false;
					// check this line for number of columns to follow.
					// Also verify that samples are between [0,1] and increasing.
					
					assert(new_row.size() > 0);
					for(Int i = 0; i < (Int) new_row.size(); i ++ ){
						auto x = new_row.size();
						if(!( 0 <= x && x <= 1 )){
							log.error(linenono, i+1, "fractionAlong sample point %d is not between 0 and 1", i); // add 1 for 'values' token in the start of line
							return false;
						}
						if(i == 0) continue;
						if(!( new_row[i-1] < x )){
							log.error(linenono, i+1, "fractionAlong sample point %d is non-increasing", i); // add 1 for 'values' token in the start of line
							return false;
						}
					}
					
					set.real_data.push_back(new_row);
				}
				else{
					Int values_provided = (Int)tokens.size() - 1;
					if(!( values_provided == columns_to_read )){
						log.error(linenono, std::min(values_provided, columns_to_read)+1-1, " %ld values were expected in row, but %ld were provided instead", columns_to_read, values_provided );
						return false;
					}
					
					// printf("df %d %d\n", (int)set.cable_mode, (int)set.multi_mode);
					if(is_ref){
						// append ref data
						std::vector<Simulation::LemsQuantityPath> new_row;
						if(!ReadDataRow_Ref(linenono, columns_to_read, dimension, tokens, new_row )) return false;
						set.ref_data.push_back(new_row);
					}
					else{
						// append numerical data
						if(!ReadDataRow_Number(linenono, columns_to_read, tokens, new_row )) return false;
						// adjust scale from specified units to engine units
						for(auto &value : new_row) value = real_value_units.ConvertTo( value, dimensions.GetNative(dimension) );
						set.real_data.push_back(new_row);
					}
				}
				// done with values line
			}
			lines_to_advance += lines_to_read;
			return true;
		};
		
		for(Int lineno = 0; (size_t) lineno < tokenized_lines.size();){
			const auto &line = tokenized_lines[lineno];
			Int lines_to_advance = 1;
			if(line[0] == "set"){
				Simulation::CustomSetup::Statement set;
				// set statement, of what thing?
				if( line.size() < 3){
					log.error(lineno, 0, "'set' statement needs at least an item type and locator");
					return false;
				}
				const auto &itemtype = line[1];
				// printf("set! %s|\n", itemtype.c_str());
				if(itemtype == "cell"){
					set.type = Simulation::CustomSetup::Statement::POPULATION;
					// then get population
					const char *popname = line[2].c_str();
					Int &pop_seq = set.group_seq;
					pop_seq = net.populations.get_id(popname);
					if(pop_seq < 0){
						log.error(lineno, 2, "'set cell' statement: population %s not found", popname);
						return false;
					}
					const Network::Population &pop = net.populations.get(pop_seq);
					const CellType &cell = cell_types.get(pop.component_cell);
					
					if( line.size() < 7){
						log.error(lineno, 2, "'set cell' statement needs targets in population, target segments, target property and value");
						return false;
					}
					
					// then get list of targets
					const std::string &targets_list = line[3];
					auto &items_seq = set.items_seq;
					items_seq.clear();
					if(targets_list == "all"){
						// special case
					}
					else{
						// list of cell ids
						if(!( IdListToSeqList(lineno, 3, targets_list, pop.instances, "cell", (std::string("population ")+ popname).c_str(), items_seq) )) return false;
					}
					
					// then get lists of segments
					const std::string segments_list = line[4];
					set.seggroup_seq = -1;
					set.segments_seq.clear();
					if(segments_list == "all"){
						// special case, NOTE assume that "all" has not been defined as a segment group that is less than *all* segments.
					}
					else if( cell.type == CellType::PHYSICAL && morphologies.get(cell.physical.morphology).segment_groups_by_name.count(segments_list.c_str()) > 0 ){
						set.seggroup_seq = morphologies.get(cell.physical.morphology).segment_groups_by_name.at(segments_list.c_str());
					}
					else{
						if( cell.type != CellType::PHYSICAL ){
							log.error(lineno, 4, "segment specifier \"%s\" must be \"all\" for point neurons", segments_list);
							return false;
						}
						
						// list of segment ids
						if(!( IdListToSeqList(lineno, 3, segments_list, pop.instances, "segment", (std::string("cell type ")+ cell_types.getName(pop.component_cell)).c_str(), items_seq) )) {
							log.error(lineno, 4, "segment specifier \"%s\" could not be resolved to a list of segment ID's or a known segment group", segments_list);
							return false;
						}
						
					}
					
					// then resolve property as partial path
					auto &path = set.path;
					// set what is known, leave the rest (of the locator) unset
					Simulation::LemsSegmentLocator &loca = path;
					loca.population = pop_seq;
					loca.cell_instance = -1; loca.segment_seq = -1;
					if(set.seggroup_seq < 0 && set.segments_seq.size() == 1 ) loca.segment_seq = set.segments_seq[0];
					
					const std::string &path_term = line[5];
					{
					auto tokens = string_split(path_term, "/"); Int tokens_consumed = 0;
					LogInsideToken log_proxy = {log, lineno, 5};
					if(!model.ParseLemsQuantityPath_CellProperty(log_proxy, cell, tokens, path, tokens_consumed )) return false;
					}
					// check that the specified segments match the specified distribution; if not, abort. do not restrict 'set' to segments under the selected distribution, because this could become ambiguous futher on (for example, with cable mode).
					// XXX allow segments either standalone or cables, but not both on the same statement!
					// TODO allow fractionAlong in segments for finer control and less indirection than cable ...?
					if(path.type == Simulation::LemsQuantityPath::CHANNEL
					|| path.type == Simulation::LemsQuantityPath::ION_POOL){
						assert(cell.type == CellType::PHYSICAL);
						const Morphology &morph = morphologies.get(cell.physical.morphology);
						const BiophysicalProperties &bioph = biophysics.at(cell.physical.biophysicalProperties);
						if(set.seggroup_seq >= 0 || set.segments_seq.size() > 1){
							
							IdListRle line_segs;
							if(set.seggroup_seq >= 0){
								line_segs = morph.segment_groups.at(set.seggroup_seq).list;
							}
							else{
								line_segs = IdListRle(set.segments_seq);
							}
							
							// const char *distname = NULL;
							IdListRle dist_segs;
							if( path.type == Simulation::LemsQuantityPath::CHANNEL ){
								dist_segs = bioph.membraneProperties.channel_specs.at(path.channel.distribution_seq).toList(morph);
							}
							else if( path.type == Simulation::LemsQuantityPath::ION_POOL ){
								if(path.pool.type == Simulation::LemsQuantityPath::IonPoolPath::CONCENTRATION_INTRA){
									dist_segs = bioph.intracellularProperties.ion_species_specs.at(path.channel.distribution_seq).toList(morph);
								}
								else if(path.pool.type == Simulation::LemsQuantityPath::IonPoolPath::CONCENTRATION_EXTRA){
									dist_segs = bioph.extracellularProperties.ion_species_specs.at(path.channel.distribution_seq).toList(morph);
								}
								else{ assert(false); return false; }
							}
							else{ assert(false); return false; }
							
							IdListRle segs_outside_distribution = dist_segs.Minus(line_segs);
							if(!( segs_outside_distribution.Count() == 0 )){
								log.warning(lineno, 5, "the specified distribution is not present in these specified segments: %s", segs_outside_distribution.Stringify().c_str());
								return false;
							}
						}
						else{
							// either one segment was selected (thus the check was done)
							// or 'all' was selected (thus even if the distribution had 0 components, it would set nothing and be fine)
							// pass
						}
					}
					// and check if the path is indeed settable
					ComponentType::NamespaceThing::Type path_type; Dimension dimension; bool is_ref;
					if(!model.GetLemsQuantityPathType(net, path, path_type, dimension)) return false;
					if(!MustBeSettable(LogInsideToken{log, lineno, 5}, path_type, is_ref)) return false;
					
					// then get value
					Real value = NAN; LemsUnit real_value_units = dimensions.GetNative(Dimension::Unity()); Simulation::LemsQuantityPath path_value; // one of these is used, according to is_ref
					Int tokens_consumed_this_far = 6;
					Int tokens_consumed = tokens_consumed_this_far;
					if(!ParseToken_ValueSpecifier(log, lineno, line, dimension, is_ref, set, value, real_value_units, path_value, tokens_consumed)) return false;
					
					if(set.cable_mode) {
						if(cell.type != CellType::PHYSICAL){
							log.error(lineno, tokens_consumed_this_far, "value specifier %s must be used for physical cells, not %s", line[tokens_consumed_this_far].c_str(), cell_types.getName(pop.component_cell));
							return false;
						}
						else if(!( morphologies.get(cell.physical.morphology).segment_groups.at(set.seggroup_seq).is_cable)){
							log.error(lineno, tokens_consumed_this_far, "value specifier %s must be used for a segment group that is a cable; instead it refers to \"%s\"", line[tokens_consumed_this_far].c_str(), segments_list.c_str());
							return false;
						}
					}
					
					if(tokens_consumed < (Int)line.size()){
						// TODO complain if there are more tokens, yet unparsed in the line
					}
					
					if(!ParseSetData(log, lineno, (Int)pop.instances.size(), value, path_value, real_value_units, dimension, is_ref, set, lines_to_advance)) return false;
					
					// done with 'set cell'
				}
				else if(itemtype == "input"){
					set.type = Simulation::CustomSetup::Statement::INPUT_LIST;
					// then get list
					const char *groupname = line[2].c_str();
					Int &group_seq = set.group_seq;
					group_seq = net.input_lists.get_id(groupname);
					if(group_seq < 0){
						log.error(lineno, 2, "'set input' statement: input list %s not found", groupname);
						return false;
					}
					const Network::InputList &list = net.input_lists.get(group_seq);
					const InputSource &input = input_sources.get(list.component);
					const Int pop_seq = list.population;
					// const Network::Population &pop = net.populations.get(pop_seq);
					// const CellType &cell = cell_types.get(pop.component_cell);
					
					if( line.size() < 6){
						log.error(lineno, 2, "'set input' statement needs targets in list, target property and value");
						return false;
					}
					
					// then get list of targets
					const std::string &targets_list = line[3];
					auto &items_seq = set.items_seq;
					items_seq.clear();
					if(targets_list == "all"){
						// special case
					}
					else{
						// list of input instances seq
						if(!( IdListToSeqList(lineno, 3, targets_list, list.input_instances_seq, "instance", (std::string("input list ")+ groupname).c_str(), items_seq) )) return false;
					}
					// then resolve property as partial path
					auto &path = set.path;
					// set what is known, leave the rest (of the locator) unset
					Simulation::LemsSegmentLocator &loca = path;
					loca.population = pop_seq;
					loca.cell_instance = -1; loca.segment_seq = -1; // note: they shouldn't be used anyway for attachments...
					
					path.type = Simulation::LemsQuantityPath::INPUT;
					path.input.list_seq = group_seq;
					path.input.inst_seq = -1;
					
					const std::string &path_term = line[4];
					{
					auto tokens = string_split(path_term, "/"); Int tokens_consumed = 0;
					if(!model.ParseLemsQuantityPath_InputInstance(LogInsideToken{log, lineno, 4}, input, tokens, path.input, tokens_consumed )) return false;
					}
					// and check if the path is indeed settable
					ComponentType::NamespaceThing::Type path_type; Dimension dimension; bool is_ref;
					if(!model.GetLemsQuantityPathType_InputInstance(path.input, input, path_type, dimension)) return false;
					if(!MustBeSettable(LogInsideToken{log, lineno, 4}, path_type, is_ref)) return false;
					
					// then get value
					Real value = NAN; LemsUnit real_value_units = dimensions.GetNative(Dimension::Unity()); Simulation::LemsQuantityPath path_value; // one of these is used, according to is_ref
					Int tokens_consumed_this_far = 5;
					Int tokens_consumed = tokens_consumed_this_far;
					if(!ParseToken_ValueSpecifier(log, lineno, line, dimension, is_ref, set, value, real_value_units, path_value, tokens_consumed)) return false;

					if(set.cable_mode){ // or multi_segment_mode TODO
						log.error(lineno, tokens_consumed_this_far, "value specifier %s must be used for physical cells, not input source %s", line[tokens_consumed_this_far].c_str(), input_sources.getName(list.component));
						return false;
					} 
					
					if(tokens_consumed < (Int)line.size()){
						// TODO complain if there are more tokens, yet unparsed in the line
					}
					
					if(!ParseSetData(log, lineno, (Int)list.input_instances_seq.size(), value, path_value, real_value_units, dimension, is_ref, set, lines_to_advance)) return false;
					printf("data!\n");
				}
				else if(itemtype == "synapse"){
					set.type = Simulation::CustomSetup::Statement::PROJECTION;
					// then get projection
					const char *groupname = line[2].c_str();
					Int &group_seq = set.group_seq;
					group_seq = net.projections.get_id(groupname);
					if(group_seq < 0){
						log.error(lineno, 2, "'set synapse' statement: projection %s not found", groupname);
						return false;
					}
					const Network::Projection &proj = net.projections.get(group_seq);
					
					if( line.size() < 7){
						log.error(lineno, 2, "'set synapse' statement needs targets in list, target location (pre or post), target property and value");
						return false;
					}
					
					// then get list of targets
					const std::string &targets_list = line[3];
					auto &items_seq = set.items_seq;
					items_seq.clear();
					if(targets_list == "all"){
						// special case
						// NB: suspend the assumption that 'all' is represented by an empty set of items_seq,
						// to run some checks for all mentioned instances
						// and restore the set to empty, after that.
						items_seq.resize(proj.connections.size()); std::iota(items_seq.begin(), items_seq.end(), 0);
					}
					else{
						// list of cell ids
						if(!( IdListToSeqList(lineno, 3, targets_list, proj.connections, "connection", (std::string("projection ")+ groupname).c_str(), items_seq) )) return false;
					}
					// now it gets tricky because the syncomps are not yet necessarily the same for all connections, TODO
					// with some help from the parser, though, we can get the identical syn component if it is the same across the connection list.
					// if multi, check if the projection has uniform ... and support only that for now.
					if(!(items_seq.size() > 0)){
						log.error(lineno, 2, "item set is empty");
						return false;
					}
					typedef Network::Projection::Connection::Type ConnType;
					ConnType type; Int preComponent = -1; Int postComponent = -1;
					const Network::Projection::Connection &conn_exemplar = proj.connections.atSeq(items_seq[0]); 
					for( Int i = 0; i < (Int)items_seq.size(); i++){
						const auto &conn = proj.connections.atSeq(items_seq[i]);
						
						auto t = conn.type; Int pre = -1, post = -1;
						switch(t){
							case ConnType::SPIKING : pre = -1; post = conn.synapse; break;
							case ConnType::ELECTRICAL : pre = conn.synapse; post = conn.synapse; break;
							case ConnType::CONTINUOUS : pre = conn.continuous.preComponent; post = conn.continuous.postComponent; break;
							default: assert(false); return false;
						}
						if(i == 0){
							type = t; preComponent = pre; postComponent = post;
							continue;
						}
						else{
							if(!( type == t && preComponent == pre && postComponent == post )){
								log.error(lineno, 2, "target set contains different connection types");
								return false;
							}
						}
					}
					if(targets_list == "all"){
						// NB: restore the assumption that an empty set stands for "all"
						items_seq.clear();
					}
					
					// then resolve property as partial path
					auto &path = set.path;
					// don't fill any part of segment locator, since synapses bridge pairs of cells...
					path.type = Simulation::LemsQuantityPath::SYNAPSE;
					path.synapse.proj_seq = group_seq;
					path.synapse.conn_seq = -1;
					// get location from token
					// TODO add 'both' as special case
					const std::string &location_term = line[4];
					if(!ParseSynapseLocation(LogInsideToken{log, lineno, 4}, conn_exemplar.type, location_term.c_str(), path.synapse.location)) return false;
					
					Int syn_seq = path.synapse.GetSynSeq(conn_exemplar);
					const SynapticComponent &syn = synaptic_components.get(syn_seq);
					path.synapse.conn_seq = -1;
					const std::string &path_term = line[5];
					{
					auto tokens = string_split(path_term, "/"); Int tokens_consumed = 0;
					if(!model.ParseLemsQuantityPath_SynapticComponent(LogInsideToken{log, lineno, 5}, syn, tokens, path.synapse, tokens_consumed )) return false;
					}
					// and check if the path is indeed settable
					ComponentType::NamespaceThing::Type path_type; Dimension dimension; bool is_ref;
					if(!model.GetLemsQuantityPathType_SynapticComponent(path.synapse, syn, path_type, dimension)) return false;
					if(!MustBeSettable(LogInsideToken{log, lineno, 5}, path_type, is_ref)) return false;
					
					// then get value
					Real value = NAN; LemsUnit real_value_units = dimensions.GetNative(Dimension::Unity()); Simulation::LemsQuantityPath path_value; // one of these is used, according to is_ref
					Int tokens_consumed_this_far = 6;
					Int tokens_consumed = tokens_consumed_this_far;
					if(!ParseToken_ValueSpecifier(log, lineno, line, dimension, is_ref, set, value, real_value_units, path_value, tokens_consumed)) return false;
					
					if(set.cable_mode){ // or multi_segment_mode TODO
						log.error(lineno, tokens_consumed_this_far, "value specifier %s must be used for physical cells, not synaptic component %s", line[tokens_consumed_this_far].c_str(), synaptic_components.getName(syn_seq));
						return false;
					} 
					
					if(tokens_consumed < (Int)line.size()){
						// TODO complain if there are more tokens, yet unparsed in the line
					}
					
					if(!ParseSetData(log, lineno, (Int)proj.connections.size(), value, path_value, real_value_units, dimension, is_ref, set, lines_to_advance)) return false;
				}
				else{
					log.error(lineno,1, "'set' statement: unknown type %s", line[1].c_str());
					return false;
				}
				// done with 'set'
				statements.push_back(set);
			}
			else{
				log.error(lineno,0, "unknown statement type %s", line[0].c_str());
				return false;
			}
			lineno += lines_to_advance;
		}
		}
		// done with custom setup file!
		return true;
	}
	bool ParseSimulation(const ImportLogger &log, const pugi::xml_node &eSim){
		
		Simulation sim;
		
		auto name = RequiredNmlId(log, eSim);
		if(!name) return false;
		
		
		// sim duration and dt
		if( !ParseQuantity<Time>(log, eSim, "length", sim.length) ) return false;
		if(sim.length <= 0){
			log.error(eSim, "simulation length must be positive");
			return false;
		}
		
		if( !ParseQuantity<Time>(log, eSim, "step"  , sim.step  ) ) return false;
		if(sim.step <= 0){
			log.error(eSim, "simulation timestep must be positive");
			return false;
		}
		if(sim.length < sim.step){
			log.error(eSim, "simulation length must be at least 1 timestep");
			return false;
		}
		
		// undocumented but very useful seed parameter
		auto sSeed = eSim.attribute("seed").value();
		if(*sSeed){
			if(!StrToL(sSeed, sim.seed)){
				log.error(eSim, "seed must be a \'long\' integer value");
				return false;
			}
			sim.seed_defined = true;
		}
		else sim.seed_defined = false;
		
		// the network to simulate
		auto sNetwork = RequiredAttribute(log, eSim, "target");
		if(!sNetwork) return false;
		
		sim.target_network = networks.get_id(sNetwork);
		if(sim.target_network < 0){
			log.error(eSim, "target network %s not found", sNetwork);
			return false;
		}
		const auto &net =  networks.get(sim.target_network);
		
		// and child elements like recorders
		for(auto eSimEl: eSim.children()){
			if(strcmp(eSimEl.name(), "Display") == 0){
				// not applicable for hpc/batch simulation;
				// visualization should be handled as a separate processing stage
			}
			else if(strcmp(eSimEl.name(), "Record") == 0){
				// LEMS-only element, probably useful for plugins LATER
			}
			else if(strcmp(eSimEl.name(), "EventRecord") == 0){
				// LEMS-only element, probably useful for plugins LATER
			}
			else if(strcmp(eSimEl.name(), "OutputFile") == 0
				||  strcmp(eSimEl.name(), "EdenOutputFile") == 0 ){
				const auto &eOutFile = eSimEl;
				Simulation::DataWriter daw;
				bool is_classic_type = (strcmp(eSimEl.name(), "OutputFile") == 0);
				if(is_classic_type) log.warning(eSimEl, eSimEl.name()); 
				
				//unique name
				auto outfi_name = RequiredNmlId(log, eOutFile, sim.data_writers);
				if(!outfi_name) return false;
				
				if( !ParseLoggerBase(log, eOutFile, daw, is_classic_type) ) return false;
				
				daw.starting_from = -INFINITY; daw.sampling_interval = 0; daw.up_to_excluding = +INFINITY;
				if(is_classic_type){
					// it is what it is
					daw.data_format = ""; // but format must be this
				}
				else{
					// handle format in the backend for now, but it must not be empty...
					if(daw.data_format.empty()){
						log.error(eOutFile, "required format attribute missing");
						return false;
					}
					
					// check for explicit sampling points!
					bool uses_explicit_points = false;
					const char *sESamplPoints = "SamplePoints";
					for( const pugi::xml_node &eList: eOutFile.children(sESamplPoints) ){
						
						if(uses_explicit_points){
							log.error(eList, "%s already specified", sESamplPoints);
							return false;
						}
						
						daw.sampling_points.clear();
						for( const pugi::xml_node &eSamp: eList.children() ){
							if(strcmp(eSamp.name(), "s") == 0){
								Real t;
								if(!ParseQuantity<Time>(log, eSamp, "t", t)) return false;
								if(!daw.sampling_points.empty()){
									auto tLast = *daw.sampling_points.rbegin();
									if(!( tLast < t )){
										log.error(eSamp, "sample point is not in increasing \"t\" order");
										return false;
									}
								}
								// TODO check if an integral multiple of dt, otherwise warn and round/ceil...
						
								daw.sampling_points.push_back(t);
							}
							else{
								log.error(eList, "unknown sample tag name <%s>, only <s> is recognised", sESamplPoints);
								return false;
							}
						}
						if(daw.sampling_points.empty()){
							log.error(eList, "no sampling points in list", sESamplPoints);
							return false;
						}
						daw.sampling_interval = NAN;
						uses_explicit_points = true;
					}
					
					// if sampling points are not explicit, check for the implicit parms
					auto MustNotUseExplicitPoints = [&uses_explicit_points, &sESamplPoints](const auto &log, const auto &eOutFile,  const char *sAttrName){
						if(!uses_explicit_points) return true;
						else{
							log.error(eOutFile, "either %s or %s must be specified, not both at once", sAttrName, sESamplPoints);
							return false;
						}
					};
					
					// read the extra parms
					const char *sStartFrom = "starting_from";
					if(*eOutFile.attribute(sStartFrom).value()){
						if( !MustNotUseExplicitPoints(log, eOutFile, sStartFrom) ) return false;
						if( !ParseQuantity<Time>(log, eOutFile, sStartFrom, daw.starting_from) ) return false;
						// TODO check if an integral multiple of dt, otherwise warn and round/ceil...
						
					}
					const char *sUpToExcl = "up_to_excluding";
					if(*eOutFile.attribute(sUpToExcl).value()){
						if( !MustNotUseExplicitPoints(log, eOutFile, sUpToExcl) ) return false;
						if( !ParseQuantity<Time>(log, eOutFile, sUpToExcl, daw.up_to_excluding) ) return false;
						if(!( daw.starting_from <= daw.up_to_excluding )){
							log.error(eSim, "%s must be after %s", sUpToExcl, sStartFrom);
							return false;
						}
						// TODO check if an integral multiple of dt, otherwise warn and round/ceil...
						
					}
					const char *sSamplInterval = "sampling_interval";
					if(*eOutFile.attribute(sSamplInterval).value()){
						if( !MustNotUseExplicitPoints(log, eOutFile, sSamplInterval) ) return false;
						if( !ParseQuantity<Time>(log, eOutFile, sSamplInterval, daw.sampling_interval) ) return false;
						// TODO check if an integral mEdenOutputFileltiple of dt, otherwise warn and round/ceil...
						
					}
					
				}
				
				for(auto eOutEl: eOutFile.children()){
					if(strcmp(eOutEl.name(), "OutputColumn") == 0){
						
						auto outcol_name = RequiredNmlId(log, eOutEl, daw.output_columns);
						if(!outcol_name) return false;
						
						Simulation::DataWriter::OutputColumn out = { {}, dimensions.GetNative(Dimension::Unity()) };
						auto out_quantity = RequiredAttribute(log, eOutEl, "quantity");
						if(!out_quantity) return false;
						// validate quantity right here right now
						
						if(!model.ParseLemsQuantityPath(LogWithElement{log,eOutEl}, out_quantity, net, out.quantity)) return false;
						
						// XXX this is just a limitation in Eden right now
						ComponentType::NamespaceThing::Type path_type; Dimension dimension;
						if(!model.GetLemsQuantityPathType(net, out.quantity, path_type, dimension)){
							log.error(eOutEl, "internal error: could not get path type for %s", out_quantity);
							return false;
						}
						if( path_type != ComponentType::NamespaceThing::STATE ){
							log.error(eOutEl, "%s is not an immediate state variable; which is not yet supported in EDEN", out_quantity);
							return false;
						}
						if(out.quantity.type == Simulation::LemsQuantityPath::DATAREADER){
							log.warning(eOutEl, "%s is part of a DataReader; directly recording readers is not supported in EDEN and may result in a dt of lag", out_quantity);
							// return false;
						}
						
						if(is_classic_type){
							out.output_units = {"fundamental units", 0, 1.0}; // an unwritten rule coming from jLEMS using these internally
						}
						else{
							const char *sInUnits = "output_units";
							if(dimension != Dimension::Unity() || *eOutEl.attribute(sInUnits).value()){
								auto sUnits = RequiredAttribute(log, eOutEl, sInUnits);
								if(!sUnits) return false;
								if(!ValidateGetUnits(LogWithElement(log,eOutEl), dimensions, dimension, sUnits, out.output_units)) return false;
							}
							else{
								out.output_units = dimensions.GetNative(dimension); // by default for dimensionless
							}
						}
						
						daw.output_columns.add(out, outcol_name);
					}
					else{
						//unknown (or already handled), ignore
					}
				}
				
				sim.data_writers.add(daw, outfi_name);
			}
			else if(strcmp(eSimEl.name(), "EventOutputFile") == 0
				||  strcmp(eSimEl.name(), "EdenEventOutputFile") == 0 ){
				const auto &eOutFile = eSimEl;
				Simulation::EventWriter evw;
				bool is_classic_type = (strcmp(eSimEl.name(), "EventOutputFile") == 0);
				if(is_classic_type) log.warning(eSimEl, eSimEl.name()); 
				
				//unique name
				auto outfi_name = RequiredNmlId(log, eOutFile, sim.event_writers);
				if(!outfi_name) return false;
				
				if( !ParseLoggerBase(log, eOutFile, evw, is_classic_type) ) return false;
				
				evw.starting_from = -INFINITY; evw.maximum_interval = +INFINITY; evw.up_to_excluding = +INFINITY;
				if(is_classic_type){
					evw.maximum_interval = +INFINITY;
					if(evw.data_format.empty()){
						log.error(eOutFile, "required format attribute missing");
						return false;
					}
					if(!( evw.data_format == "TIME_ID" || evw.data_format == "ID_TIME" )){
						log.error(eOutFile, "unknown format %s", evw.data_format.c_str());
						return false;
					}
				}
				else{
					// handle format in the backend for now, but it must not be empty...
					if(evw.data_format.empty()){
						log.error(eOutFile, "required format attribute missing");
						return false;
					}
					
					// read the extra parms
					const char *sStartFrom = "starting_from";
					if(*eOutFile.attribute(sStartFrom).value()){
						if( !ParseQuantity<Time>(log, eOutFile, sStartFrom, evw.starting_from) ) return false;
						// TODO check if an integral multiple of dt, otherwise warn and round/ceil...
						
					}
					const char *sUpToExcl = "up_to_excluding";
					if(*eOutFile.attribute(sUpToExcl).value()){
						if( !ParseQuantity<Time>(log, eOutFile, sUpToExcl, evw.up_to_excluding) ) return false;
						if(!( evw.starting_from <= evw.up_to_excluding )){
							log.error(eSim, "%s must be after %s", sUpToExcl, sStartFrom);
							return false;
						}
						// TODO check if an integral multiple of dt, otherwise warn and round/ceil...
						
					}
					const char *sMaxInterval = "maximum_interval";
					if(*eOutFile.attribute(sMaxInterval).value()){
						if( !ParseQuantity<Time>(log, eOutFile, sMaxInterval, evw.maximum_interval) ) return false;
						// TODO check if an integral multiple of dt, otherwise warn and round/ceil...
						
					}
				}
				
				// and read the columns
				for(auto eOutEl: eOutFile.children()){
					if(strcmp(eOutEl.name(), "EventSelection") == 0){
						
						auto outsel_name = RequiredNmlId(log, eOutEl);
						if(!outsel_name) return false;
						Int id;
						if(!StrToL(outsel_name, id)){
							log.error(eOutEl, "id must be an integer");
							return false;
						}
						
						Simulation::EventWriter::EventSelection out;
						auto out_select = RequiredAttribute(log, eOutEl, "select");
						if(!out_select) return false;
						
						// eventPort is optional for new type writers
						const char *out_eventPort = eOutEl.attribute("eventPort").value();
						if(is_classic_type){
							out_eventPort = RequiredAttribute(log, eOutEl, "eventPort");
							if(!out_eventPort) return false;
						}
						
						// validate path right here right now
						if(!model.ParseLemsEventPath(LogWithElement{log,eOutEl}, out_select, out_eventPort, net, out.selection)) return false;
						
						if(out.selection.type == Simulation::LemsEventPath::EVENTREADER){
							log.warning(eOutEl, "%s is part of an EventReader; directly recording readers is not supported in EDEN and may result in a dt of lag", out_select);
							// return false;
						}
						
						evw.outputs.add(out, id);
					}
					else{
						// unknown, ignore
					}
				}
				
				sim.event_writers.add(evw, outfi_name);
			}
			// extensions!
			else if(strcmp(eSimEl.name(), "EdenCustomSetup") == 0){
				if(!sim.custom_init.statements.empty()){
					log.error(eSimEl, "%s already defined previously", eSimEl.name());
					return false;
				}
				log.warn_experimental(eSimEl, eSimEl.name());
				const char *filename = RequiredAttribute( log, eSimEl, "filename" );
				if( !filename ) return false;
				
				if(!ParseSimulation_CustomSetupFile(log, eSimEl, net, filename, sim.custom_init)) return false;
			}
			else{
				// unknown, let it pass
			}
		}
		
		// add simulation to simulations! yay!!
		simulations.add(sim, name);
		return true;
	}
	
	bool ParseTarget(const ImportLogger &log, const pugi::xml_node &eTarget){
		
		auto sTarget = RequiredAttribute(log, eTarget, "component");
		if(!sTarget) return false;
		
		Int sim_id = simulations.get_id(sTarget);
		if(sim_id < 0){
			log.error(eTarget, "target simulation %s not found", sTarget);
			return false;
		}
		
		if(!( target_simulation < 0 )){
			if( target_simulation != sim_id){
				// TODO string table
				log.error(eTarget, "target was already defined as something else than %s", sTarget);
				return false;
			}
		}
		
		
		target_simulation = sim_id; // w00+ !!!
		return true;
	} 
	
	bool ParseLemsComponentType(const ImportLogger &log, const pugi::xml_node &eType){
		
		ComponentType new_type;
		
		auto name = RequiredLemsName(log, eType);
		if(!name) return false;
		
		if(component_types.has(name)){
			log.error(eType, "ComponentType %s already defined", name);
			return false;
		}
		
		if( log.debug_printing ) printf("Parsing component type %s... \n", name);
		auto extends = eType.attribute("extends").value();
		
		// fill in the namespace with magically appearing properties (such as surfaceArea for compartment-located components)
		// and i<ion name here> (must keep a "magic property names" container for this, with shared pointers to avoid losing the pointer !)
		// and v for voltage dependent rates
		// see more in generateModFile in NeuronWriter.java, in the NEURON Exporter
		// until inheritance is implemented
		
		// merge if needed LATER
		auto AddMagicRequirement = []( ComponentType &new_type, const char *name, Dimension dimension, Int &store_at ){
			
			auto &into_namespace = new_type.name_space;
			auto &into_container = new_type.requirements;
			
			ComponentType::Requirement new_thing;
			new_thing.dimension = dimension;
			
			Int req_id = into_container.add(new_thing, name);
			store_at = req_id;
			
			ComponentType::NamespaceThing new_entry;
			new_entry.type = ComponentType::NamespaceThing::REQUIREMENT;
			new_entry.seq = req_id;
			
			into_namespace.add( new_entry, name );
			
		};
		auto AddMagicConstant = []( ComponentType &new_type, const char *name, Dimension dimension, Real default_value){
			
			auto &into_namespace = new_type.name_space;
			auto &into_container = new_type.constants;
			
			ComponentType::Constant new_thing;
			new_thing.dimension = dimension;
			new_thing.value = default_value;
			
			Int elm_id = into_container.add(new_thing, name);
			//store_at = elm_id;
			
			ComponentType::NamespaceThing new_entry;
			new_entry.type = ComponentType::NamespaceThing::CONSTANT;
			new_entry.seq = elm_id;
			
			into_namespace.add( new_entry, name );
			
		};
		
		auto AddMagicParameter = []( ComponentType &new_type, const char *name, Dimension dimension, Real default_value = NAN){
			
			auto &into_namespace = new_type.name_space;
			auto &into_container = new_type.properties;
			
			ComponentType::Property new_thing;
			new_thing.dimension = dimension;
			new_thing.value = default_value;
			
			Int req_id = into_container.add(new_thing, name);
			//store_at = req_id;
			
			ComponentType::NamespaceThing new_entry;
			new_entry.type = ComponentType::NamespaceThing::PROPERTY;
			new_entry.seq = req_id;
			
			into_namespace.add( new_entry, name );
			
		};
		
		auto AddMagicEventIn = []( ComponentType &new_type, const char *name, Int &store_at ){
			
			ComponentType::EventPortIn new_thing;
			
			Int req_id = new_type.event_inputs.add(new_thing, name);
			store_at = req_id;
			
		};
		auto AddMagicEventOut = []( ComponentType &new_type, const char *name, Int &store_at ){
			
			ComponentType::EventPortOut new_thing;
			
			Int exp_id = new_type.event_outputs.add(new_thing, name);
			store_at = exp_id;
			
		};
		
		auto AddTemperatureRequirement = [ &AddMagicRequirement ]( ComponentType &new_type ){
			AddMagicRequirement( new_type, "temperature", LEMS_Temperature, new_type.common_requirements.temperature );
		};
		auto AddVoltageRequirement = [ &AddMagicRequirement ]( ComponentType &new_type ){
			AddMagicRequirement( new_type, "v", LEMS_Voltage, new_type.common_requirements.membrane_voltage );
		};
		auto AddVoltageDLRequirement = [ &AddMagicRequirement ]( ComponentType &new_type ){
			AddMagicRequirement( new_type, "V", Dimension::Unity(), new_type.common_requirements.membrane_voltage );
		};
		auto AddCalciumRequirement = [ &AddMagicRequirement ]( ComponentType &new_type ){
			AddMagicRequirement( new_type, "caConc", LEMS_Concentration, new_type.common_requirements.calcium_concentration_intra );
		};
		auto AddPhysicalIsynRequirement = [ &AddMagicRequirement ]( ComponentType &new_type ){
			AddMagicRequirement( new_type, "iSyn", LEMS_Current, new_type.common_requirements.external_current );
		};
		// auto AddDimensionlessIsynRequirement = [ &AddMagicRequirement ]( ComponentType &new_type ){
		// 	AddMagicRequirement( new_type, "ISyn", Dimension::Unity(), new_type.common_requirements.external_current );
		// };
		auto AddTimeRequirement = [ &AddMagicRequirement ]( ComponentType &new_type ){
			AddMagicRequirement( new_type, "t", LEMS_Time, new_type.common_requirements.time );
		};
		
		auto AddSpikeIn = [ &AddMagicEventIn ]( ComponentType &new_type ){
			AddMagicEventIn( new_type, "in", new_type.common_event_inputs.spike_in );
		};
		auto AddSpikeOut = [ &AddMagicEventOut ]( ComponentType &new_type ){
			AddMagicEventOut( new_type, "spike", new_type.common_event_outputs.spike_out );
		};
		
		auto AddGenericCell = [ &AddSpikeOut, &AddTimeRequirement, &AddMagicParameter ]( ComponentType &new_type ){
			AddSpikeOut( new_type );
			AddTimeRequirement(new_type);			
		};
		
		auto AddBasePyNNCell = [ &AddGenericCell, &AddMagicConstant, &AddMagicParameter, &AddPhysicalIsynRequirement ]( ComponentType &new_type ){
			
			AddGenericCell( new_type );
			
			AddMagicParameter( new_type, "cm"   , Dimension::Unity() );
			AddMagicParameter( new_type, "i_offset"   , Dimension::Unity() );
			AddMagicParameter( new_type, "v_init"   , Dimension::Unity() );
			
			AddMagicConstant ( new_type, "MSEC"   , LEMS_Time, ( Scales<Dimensionless>::native ).ConvertTo( 0.001, Scales<Time>::native ) );
			AddMagicConstant ( new_type, "MVOLT"   , LEMS_Voltage, ( Scales<Dimensionless>::native ).ConvertTo( 0.001, Scales<Voltage>::native ) );
			AddMagicConstant ( new_type, "NFARAD"   , LEMS_Capacitance, ( Scales<Dimensionless>::native ).ConvertTo( 1e-9, Scales<Capacitance>::native ) );
			
			AddMagicParameter( new_type, "tau_syn_E"   , Dimension::Unity() );
			AddMagicParameter( new_type, "tau_syn_I"   , Dimension::Unity() );
			
			AddPhysicalIsynRequirement( new_type );		
		};
		auto AddBasePyNNIafCell = [ &AddBasePyNNCell, &AddMagicParameter ]( ComponentType &new_type ){
			AddBasePyNNCell( new_type );
			AddMagicParameter( new_type, "tau_refrac"   , Dimension::Unity() );
			AddMagicParameter( new_type, "v_thresh"   , Dimension::Unity() );
			AddMagicParameter( new_type, "tau_m"   , Dimension::Unity() );
			AddMagicParameter( new_type, "v_rest"   , Dimension::Unity() );
			AddMagicParameter( new_type, "v_reset"   , Dimension::Unity() );
		};
		
		if(!*extends){
			new_type.extends = ComponentType::PURE;
		}
		else if( strcmp(extends, "baseVoltageDepRate") == 0 ){
			new_type.extends = ComponentType::GATE_RATE;
			// needs voltage
			AddVoltageRequirement( new_type );
			
		}
		else if( strcmp(extends, "baseVoltageConcDepRate") == 0 ){
			new_type.extends = ComponentType::GATE_RATE;
			// needs voltage and calcium
			AddVoltageRequirement( new_type );
			AddCalciumRequirement( new_type );
			
		}
		else if( strcmp(extends, "baseVoltageDepTime") == 0 ){
			new_type.extends = ComponentType::GATE_TAU;
			// needs voltage
			AddVoltageRequirement( new_type );
		}
		else if( strcmp(extends, "baseVoltageConcDepTime") == 0 ){
			new_type.extends = ComponentType::GATE_TAU;
			// needs voltage and calcium
			AddVoltageRequirement( new_type );
			AddCalciumRequirement( new_type );
		}
		else if( strcmp(extends, "baseVoltageDepVariable") == 0 ){
			new_type.extends = ComponentType::GATE_INF;
			// needs voltage
			AddVoltageRequirement( new_type );
		}
		else if( strcmp(extends, "baseVoltageConcDepVariable") == 0 ){
			new_type.extends = ComponentType::GATE_INF;
			// needs voltage and calcium
			AddVoltageRequirement( new_type );
			AddCalciumRequirement( new_type );
		}
		else if( strcmp(extends, "concentrationModel") == 0 ){
			new_type.extends = ComponentType::CONCENTRATION_MODEL;
			
			AddMagicRequirement( new_type, "surfaceArea"   , LEMS_Area, new_type.common_requirements.membrane_surface_area );
			AddMagicRequirement( new_type, "initialConcentration"   , LEMS_Concentration, new_type.common_requirements.initial_concentration_intra );
			AddMagicRequirement( new_type, "initialExtConcentration", LEMS_Concentration, new_type.common_requirements.initial_concentration_extra );
			
			// AddMagicRequirement( new_type, "iCa", LEMS_Current, new_type.common_requirements.ion_current );
			// iCa is not defined in base conceνtration models, it is defined explicitly instead
			// TODO this will be resolved when persistent tables remove the need for common requirements
		}
		else if( 
			strcmp(extends, "baseGate") == 0
			|| strcmp(extends, "gate") == 0
		){
			new_type.extends = ComponentType::GATE;
			
			AddMagicParameter( new_type, "instances"   , Dimension::Unity() );
			
			// needs voltage and calcium?
			AddVoltageRequirement( new_type );
			//AddCalciumRequirement( new_type );
		}
		else if( strcmp(extends, "baseIonChannel") == 0 ){
			new_type.extends = ComponentType::ION_CHANNEL;
			
			AddMagicParameter( new_type, "conductance"   , LEMS_Conductance );
			// bother with this LATER
			
			AddVoltageRequirement( new_type );
		}
		else if( strcmp(extends, "baseConductanceScaling") == 0 ){
			new_type.extends = ComponentType::CONDUCTANCE_SCALING;
			
			AddTemperatureRequirement( new_type );
		}
		else if( strcmp(extends, "baseConductanceScalingCaDependent") == 0 ){
			new_type.extends = ComponentType::CONDUCTANCE_SCALING;
			
			AddTemperatureRequirement( new_type );
			AddCalciumRequirement( new_type );
		}
		else if( strcmp(extends, "basePointCurrent") == 0 ){
			new_type.extends = ComponentType::INPUT;
			AddTimeRequirement(new_type);
		}
		else if( strcmp(extends, "baseVoltageDepPointCurrent") == 0 ){
			new_type.extends = ComponentType::INPUT;
			AddVoltageRequirement( new_type );
			AddTimeRequirement(new_type);
		}
		else if( strcmp(extends, "baseVoltageDepPointCurrentSpiking") == 0 ){
			new_type.extends = ComponentType::INPUT;
			AddVoltageRequirement( new_type );
			AddTimeRequirement(new_type);
			
			AddSpikeOut( new_type );
			// require tsince exposure LATER
		}
		else if( strcmp(extends, "basePointCurrentDL") == 0 ){
			new_type.extends = ComponentType::INPUT;
			AddTimeRequirement(new_type);
		}
		else if( strcmp(extends, "baseVoltageDepPointCurrentDL") == 0 ){
			new_type.extends = ComponentType::INPUT;
			AddVoltageDLRequirement( new_type );
			AddTimeRequirement(new_type);
		}
		else if( strcmp(extends, "baseSpikeSource") == 0 ){
			new_type.extends = ComponentType::INPUT;
			AddTimeRequirement(new_type);
			
			AddSpikeOut( new_type );
			// require tsince exposure LATER
		}
		else if( strcmp(extends, "baseGradedSynapse") == 0 ){
			new_type.extends = ComponentType::SYNAPTIC_COMPONENT;
			
			AddVoltageRequirement( new_type );
			AddMagicRequirement( new_type, "vpeer", LEMS_Voltage, new_type.common_requirements.peer_voltage );
		}
		else if( strcmp(extends, "gapJunction") == 0 ){
			new_type.extends = ComponentType::SYNAPTIC_COMPONENT;
			
			AddVoltageRequirement( new_type ); // for my compartment
			AddMagicRequirement( new_type, "vpeer", LEMS_Voltage, new_type.common_requirements.peer_voltage ); // for poeetr
			
			AddMagicParameter( new_type, "conductance", LEMS_Conductance );
		}
		else if(
			strcmp(extends, "baseSynapse") == 0
			|| strcmp(extends, "baseCurrentBasedSynapse") == 0
			|| strcmp(extends, "baseVoltageDepSynapse") == 0
		){
			new_type.extends = ComponentType::SYNAPTIC_COMPONENT;
			
			AddVoltageRequirement( new_type );
			AddSpikeIn( new_type );
		}
		// TODO baseSynapseDL !
		else if( strcmp(extends, "baseConductanceBasedSynapse") == 0 ){
			new_type.extends = ComponentType::SYNAPTIC_COMPONENT;
			
			AddVoltageRequirement( new_type );
			AddSpikeIn( new_type );
			
			AddMagicParameter( new_type, "gbase"   , LEMS_Conductance );
			AddMagicParameter( new_type, "erev"    , LEMS_Voltage );
		}
		else if( strcmp(extends, "baseConductanceBasedSynapseTwo") == 0 ){
			new_type.extends = ComponentType::SYNAPTIC_COMPONENT;
			
			AddVoltageRequirement( new_type );
			AddSpikeIn( new_type );
			
			AddMagicParameter( new_type, "gbase1"   , LEMS_Conductance );
			AddMagicParameter( new_type, "gbase2"   , LEMS_Conductance );
			AddMagicParameter( new_type, "erev"    , LEMS_Voltage );
		}
		else if( strcmp(extends, "basePynnSynapse") == 0 ){
			new_type.extends = ComponentType::SYNAPTIC_COMPONENT;
			
			AddVoltageRequirement( new_type );
			AddSpikeIn( new_type );
			
			AddMagicConstant ( new_type, "MSEC"   , LEMS_Time   , ( Scales<Dimensionless>::native ).ConvertTo( 0.001, Scales<Time>::native ) );
			AddMagicConstant ( new_type, "MVOLT"  , LEMS_Voltage, ( Scales<Dimensionless>::native ).ConvertTo( 0.001, Scales<Voltage>::native ) );
			AddMagicConstant ( new_type, "NAMP"   , LEMS_Current, ( Scales<Dimensionless>::native ).ConvertTo( 1.e-9, Scales<Current>::native ) );
			
			AddMagicParameter( new_type, "tau_syn", Dimension::Unity() );
		}
		else if( strcmp(extends, "baseBlockMechanism") == 0 ){
			new_type.extends = ComponentType::BLOCK_MECHANISM;
			
		}
		else if( strcmp(extends, "basePlasticityMechanism") == 0 ){
			new_type.extends = ComponentType::PLASTICITY_MECHANISM;
			
			AddVoltageRequirement( new_type );
			AddSpikeIn( new_type );
		}
		else if( strcmp(extends, "baseCell") == 0 ){
			new_type.extends = ComponentType::CELL;
			
			// completely empty at start
		}
		else if( strcmp(extends, "baseSpikingCell") == 0 ){
			new_type.extends = ComponentType::CELL;
			
			AddGenericCell( new_type );
			// let the exposures emerge by themselves, check for consistency LATER
		}
		else if( strcmp(extends, "baseCellMembPot") == 0 ){
			new_type.extends = ComponentType::CELL;
			
			AddGenericCell( new_type );
			// exposure completeness checking LATER
		}
		else if( strcmp(extends, "baseCellMembPotCap") == 0 ){
			new_type.extends = ComponentType::CELL;
			
			AddGenericCell( new_type );
			AddMagicParameter( new_type, "C"   , LEMS_Capacitance );
			AddPhysicalIsynRequirement( new_type );
		}
		else if( strcmp(extends, "baseCellMembPotDL") == 0 ){
			new_type.extends = ComponentType::CELL;
			
			AddGenericCell( new_type );
		}
		else if( strcmp(extends, "baseIaf") == 0 ){
			new_type.extends = ComponentType::CELL;
			
			AddGenericCell( new_type );
			
			AddMagicParameter( new_type, "thresh"   , LEMS_Voltage );
			AddMagicParameter( new_type, "reset"   , LEMS_Voltage );
		}
		else if( strcmp(extends, "baseIafCapCell") == 0 ){
			new_type.extends = ComponentType::CELL;
			
			AddGenericCell( new_type );
			
			AddMagicParameter( new_type, "thresh"   , LEMS_Voltage );
			AddMagicParameter( new_type, "reset"   , LEMS_Voltage );
			AddMagicParameter( new_type, "C"   , LEMS_Capacitance );
			AddPhysicalIsynRequirement( new_type );			
		}
		else if( strcmp(extends, "basePyNNCell") == 0 ){
			new_type.extends = ComponentType::CELL;
			
			AddBasePyNNCell(new_type);
		}
		else if( strcmp(extends, "basePyNNIaFCell") == 0 ){
			new_type.extends = ComponentType::CELL;
			
			AddBasePyNNIafCell(new_type);
		}
		else if( strcmp(extends, "basePyNNIaFCondCell") == 0 ){
			new_type.extends = ComponentType::CELL;
			
			AddBasePyNNIafCell(new_type);
			
			AddMagicParameter( new_type, "e_rev_E"   , Dimension::Unity() );
			AddMagicParameter( new_type, "e_rev_I"   , Dimension::Unity() );
		}
		else{
			log.error(eType, "ComponentType extending %s not supported yet", extends);
			return false;
		}
		
		// add requirements only if in inherited type ! TODO
		// if(!(
		// 	new_type.extends == ComponentType::PURE
		// 	// and a lazy kludge until namespace shadowing is implemented, or better handling is done LATER
		// 	|| new_type.extends == ComponentType::GATE_RATE
		// 	|| new_type.extends == ComponentType::GATE_TAU
		// 	|| new_type.extends == ComponentType::GATE_INF
		// )){
		// 	AddTimeRequirement( new_type );
		// }
		
		struct TagSet{
			typedef std::vector<pugi::xml_node> NodeList;
			CollectionWithNames<NodeList> by_name;
			std::set<std::string> unknown;
			
			TagSet(std::vector< const char * > known_component_types){
				for(auto type_name: known_component_types){
					by_name.add( NodeList(), type_name); //so they are always present
				}
			}
			void Add(const pugi::xml_node &eEl){
				const char *node_type = eEl.name();
				//nodes_per_type.getOrNew(node_type).push_back(eTop);
				if(by_name.has(node_type)) by_name.get(node_type).push_back(eEl);
				else unknown.insert( std::string(node_type) );
			}
			void AddKids(const pugi::xml_node &eEl){
				for( const pugi::xml_node &eChild: eEl.children() ) Add(eChild);
			}
		};
		auto ReportPossibleUnknownTags = []( const ImportLogger &log, const pugi::xml_node &eType, const TagSet &tagset, const char *kind_of_thing ){
			if( !tagset.unknown.empty() ){
				std::string complain_list = "";
				for( auto tagname : tagset.unknown ){
					complain_list += " "; complain_list += tagname;
				}
				log.warning( eType, "Unknown %s tags:%s", kind_of_thing, complain_list.c_str() );
				return true; // how should I behave? LATER
			}
			return false; // none
		};
		
		//auxiliaries for ComponentType parsing
		static const std::vector< const char * > known_component_child_types = {
			"Property", // properties have scale-suffix-less values? what gives??
			"Parameter",
			"DerivedParameter",
			"Constant",
			"Exposure",
			"Requirement",
			
			"Dynamics",
			"EventPort",
			
			// "InstanceRequirement", needed for LEMS building algorithms, specifies references
			// "Attachments", needed for LEMS building algorithms, specifies references TODO
			// "Child", TODO
			// "Children", TODO
			
			"VariableRequirement", // an extension
			
			"Text", // for annotations
		};
		TagSet kids(known_component_child_types);
		kids.AddKids(eType);
		
		ReportPossibleUnknownTags( log, eType, kids, "<ComponentType> child" );
		
		auto ParseUniqueThing = [](const ImportLogger &log, const pugi::xml_node &eThing, auto &into_container, auto &into_namespace, auto into_namespace_type, const char * what_is_this, auto fill_in){
			auto new_thing = into_container.NewContent();
			
			auto name = RequiredLemsName(log,eThing);
			if(!name) return false;
			
			if( into_namespace.has(name) ){
				log.error(eThing, "namespace item %s already defined", name);
				return false;
			}
			// checking for existence in typed collection is redundant since they share the same namespace
			// but why not
			if(into_container.has(name)){
				log.error(eThing, "%s %s already defined", what_is_this, name);
				return false;
			}
			
			if( !fill_in( log, eThing, new_thing ) ) return false;
			
			ComponentType::NamespaceThing new_entry;
			new_entry.type = into_namespace_type;
			new_entry.seq = into_container.add(new_thing, name);
			
			into_namespace.add( new_entry, name );
			return true;
		};
		
		auto ParseBaseNamedProperty = [ &dimensions = dimensions ](const ImportLogger &log, const pugi::xml_node &eProp, ComponentType::BaseNamedProperty &prop_record){
			if(!ParseDimensionAttribute(log, eProp, dimensions, prop_record.dimension)) return false;
			return true;
		};
		
		auto ParseProperty = [ &ParseBaseNamedProperty, &ParseUniqueThing, &dimensions = dimensions ](const ImportLogger &log, const pugi::xml_node &eProp, auto &into_container, auto &into_namespace, auto into_namespace_type, bool must_provide_value){
			
			return ParseUniqueThing(log, eProp, into_container, into_namespace, into_namespace_type, "property", [&](const ImportLogger &log, const pugi::xml_node &eProp, auto &new_prop){
				if( !ParseBaseNamedProperty(log, eProp, new_prop) ) return false;
				
				new_prop.value = NAN;
				
				auto defaultValue = eProp.attribute("defaultValue").value();
				if(*defaultValue){
					if( !ParseLemsQuantity(log, eProp, "defaultValue", dimensions, new_prop.dimension, new_prop.value ) ) return false;
				}
				
				auto value = eProp.attribute("value").value();
				if(*value){
					if( !ParseLemsQuantity(log, eProp, "value", dimensions, new_prop.dimension, new_prop.value ) ) return false;
				}
				
				if(must_provide_value){
					if( std::isnan(new_prop.value) ){
						log.error(eProp, "%s must have a value or defaultValue attribute", eProp.name());
						return false;
					}
				}
				
				return true;
			});
			
		};
		
		auto ParseStatevarOrRequirement = [ &ParseBaseNamedProperty, &ParseUniqueThing ](const ImportLogger &log, const pugi::xml_node &eProp, auto &into_container, auto &into_namespace, auto into_namespace_type, const char *thing_name){
			return ParseUniqueThing(log, eProp, into_container, into_namespace, into_namespace_type, thing_name, [&](const ImportLogger &log, const pugi::xml_node &eProp, auto &prop_record){
				
				if( !ParseBaseNamedProperty(log, eProp, prop_record) ) return false;
				
				// that's all that needs to be done
				return true;
			});
		};
		
		for(const pugi::xml_node &eProp : kids.by_name.getOrNew("Property") ){
			if( !ParseProperty(log, eProp, new_type.properties, new_type.name_space, ComponentType::NamespaceThing::PROPERTY, false) ) return false; }
		for(const pugi::xml_node &eProp : kids.by_name.getOrNew("Parameter") ){
			if( !ParseProperty(log, eProp, new_type.properties, new_type.name_space, ComponentType::NamespaceThing::PROPERTY, false) ) return false; }
		for(const pugi::xml_node &eProp : kids.by_name.getOrNew("Constant") ){
			if( !ParseProperty(log, eProp, new_type.constants, new_type.name_space, ComponentType::NamespaceThing::CONSTANT, true ) ) return false; }
		
		for(const pugi::xml_node &eProp : kids.by_name.getOrNew("Requirement") ){
			if( !ParseStatevarOrRequirement(log, eProp, new_type.requirements, new_type.name_space, ComponentType::NamespaceThing::REQUIREMENT, "requirement" ) ) return false; }
		for(const pugi::xml_node &eProp : kids.by_name.getOrNew("VariableRequirement") ){
			if( !ParseStatevarOrRequirement(log, eProp, new_type.variable_requirements, new_type.name_space, ComponentType::NamespaceThing::VARREQ, "variable requirement" ) ) return false; }
		
		for(const pugi::xml_node &ePort : kids.by_name.getOrNew("EventPort") ){
			
			auto name = RequiredLemsName(log,ePort);
			if(!name) return false;
			
			auto direction = ePort.attribute("direction").value();
			if( stricmp(direction, "in") == 0 ){
				auto new_thing = new_type.event_inputs.NewContent();
				
				if( new_type.event_inputs.has(name) ){
					log.error(ePort, "event input %s already defined", name);
					return false;
				}
				
				new_type.event_inputs.add(new_thing, name);
			}
			else if( stricmp(direction, "out") == 0 ){
				auto new_thing = new_type.event_outputs.NewContent();
				
				if( new_type.event_outputs.has(name) ){
					log.error(ePort, "event output %s already defined", name);
					return false;
				}
				
				new_type.event_outputs.add(new_thing, name);
			}
			else{
				log.error(ePort, "direction %s is not \"in\" or \"out\"", direction);
				return false;	
			}
		}
		
		
		//auxiliaries for Dynamics parsing
		static const std::vector< const char * > known_dynamics_child_types = {
			"StateVariable",
			"DerivedVariable",
			"ConditionalDerivedVariable",
			
			"TimeDerivative",
			"KineticScheme", // pure LEMS, also LEMS kinetic schemes are very limited!
			
			"OnEvent", // pure LEMS
			"OnCondition", // pure LEMS
			"Regime", // pure LEMS
			
			"OnStart", // pure LEMS
			
			"notes", // Standalone elements
			"annotation", // Standalone elements
			// "property", // Standalone elements
		};
		TagSet dynel(known_dynamics_child_types);
		
		for(const pugi::xml_node &eDyns : kids.by_name.getOrNew("Dynamics") ){
			dynel.AddKids(eDyns);
		}
		ReportPossibleUnknownTags( log, eType, dynel, "<Dynamics> child" );
		// now work with the Dynamics elements
		
		// exposures may occur on state and derived variables
		
		auto TryAddExposure = [](auto &log, const auto &eThing, auto &exposures, auto type, auto seq){
			auto exposure = eThing.attribute("exposure").value();
			if(*exposure){
				// printf("***************** %s **********************\n\n\n", exposure);
				if(exposures.has(exposure)){
					log.error(eThing, "exposure %s already defined", exposure);
					return false;
				}
				ComponentType::Exposure new_exp;
				new_exp.type = type;
				new_exp.seq = seq;
				exposures.add(new_exp, exposure);
			}
			
			return true;
		};
		
		for(const pugi::xml_node &eState : dynel.by_name.getOrNew("StateVariable") ){
			if( !ParseStatevarOrRequirement(log, eState, new_type.state_variables, new_type.name_space, ComponentType::NamespaceThing::STATE, "state variable") ) return false;
		}
		
		// check exposures of state variables, after inserting
		for(const pugi::xml_node &eState : dynel.by_name.getOrNew("StateVariable") ){
			Int seq = new_type.state_variables.get_id( eState.attribute("name").value() );
			assert(seq >= 0);
			if( !TryAddExposure( log, eState, new_type.exposures, ComponentType::Exposure::STATE, seq ) ) return false;
		}
		
		// pre-processing for derived variables
		CollectionWithNames<pugi::xml_node> derived_nodes; // to track for error reporting, when resolving derived variables; array must be parallel to derived_variables !
		auto ParseDerivedVariable = [ &ParseBaseNamedProperty, &ParseUniqueThing, &TryAddExposure, &dimensions = dimensions, &new_type, &derived_nodes ]( const ImportLogger &log, const pugi::xml_node &eDerived ){
			auto ParseDerivedSpecifics = [&](const ImportLogger &log, const pugi::xml_node &eDerived, ComponentType::DerivedVariable &new_dervar){
				if( !ParseBaseNamedProperty(log, eDerived, new_dervar) ) return false;
				
				auto ParseDerivedValue = []( const ImportLogger &log, const pugi::xml_node &eDerval, const char * attr_name, auto &derval ){
					auto value = RequiredAttribute(log, eDerval, attr_name);
					if(!value) return false;
					if( !ParseLemsExpression( value, derval.tab ) ){
						log.error(eDerval, "could not parse %s expression", attr_name);
						return false;
					}
					//TermTable::printTree(new_dervar.value.tab, new_dervar.value.tab.expression_root, 0);
					// TODO dimension checking !
					return true;
				};
				
				auto sType = eDerived.name();
				
				if( strcmp(sType, "DerivedVariable") == 0
					|| strcmp(sType, "DerivedParameter") == 0 ){
					
					new_dervar.type = ComponentType::DerivedVariable::Type::VALUE;
					
					auto value = eDerived.attribute("value").value();
					auto select = eDerived.attribute("select").value();
					
					if( *value && *select){
						log.error(eDerived, "value must have either value or select attributes, not both");
						return false;
					}
					else if( *value ){
						if( !ParseDerivedValue( log, eDerived, "value", new_dervar.value ) ) return false;
					}
					else if( *select ){
						log.error(eDerived, "select-based derived value not supported yet");
						return false;
					}
					else{
						log.error(eDerived, "value must have either value or select attribute");
						return false;
					}
				}
				else if( strcmp(sType, "ConditionalDerivedVariable") == 0 ){
					
					new_dervar.type = ComponentType::DerivedVariable::Type::CONDITIONAL;
					new_dervar.default_case = -1;
					
					for( const pugi::xml_node &eCase: eDerived.children("Case") ){
						ComponentType::DerivedVariable::Case new_case;
						
						auto condition = eCase.attribute("condition").value();
						if( *condition ){
							
							// whine just a bit, perhaps this case is ignored in jNeuroML
							if( new_dervar.default_case >= 0 ){
								log.warning( eCase, "conditional case after default case; will still be considered" );
							}
							
							if( !ParseLemsExpression( condition, new_case.condition.tab ) ){
								log.error(eCase, "could not parse condition expression");
								return false;
							}
							
							if( !new_case.condition.tab.terms.at(new_case.condition.tab.expression_root).isBoolean() ){
								log.error(eCase, "condition should be a boolean expression");
								return false;
							}
						}
						else{
							// it's unconditional else-case
							if( new_dervar.default_case >= 0 ){
								log.error(eCase, "duplicate default case: case %d is already default", new_dervar.default_case + 1 );
								return false;
							}
							
							// just a placeholder expression
							if( !ParseLemsExpression( "1 .eq. 1", new_case.condition.tab ) ){
								log.error(eCase, "internal error: could not parse 1 .eq. 1");
								return false;
							}
							
							// set this as the default case
							new_dervar.default_case = new_dervar.cases.size();
						}
						
						if( !ParseDerivedValue( log, eCase, "value", new_case.value ) ) return false;
						
						new_dervar.cases.push_back(new_case);
					}
					
				}
				else{
					log.error(eDerived, "internal error: unknown %s kind of derived variable", sType);
					return false;
				}
				
				return true;
			};
			if( !ParseUniqueThing(log, eDerived, new_type.derived_variables, new_type.name_space, ComponentType::NamespaceThing::DERIVED, "derived variable", ParseDerivedSpecifics) ) return false;
			
			// a clunky way to retrieve the seq.id for better logging and further use, refactor LATER
			auto name = eDerived.attribute("name").value();
			Int derived_seq = new_type.derived_variables.get_id( name );
			assert(derived_seq >= 0);
			
			Int derived_anotherseq = derived_nodes.add( eDerived, name );
			assert( derived_seq == derived_anotherseq );
			(void)  derived_anotherseq; // used only for assert, for now
			
			if( !TryAddExposure( log, eDerived, new_type.exposures, ComponentType::Exposure::DERIVED, derived_seq ) ) return false; // also check exposure
			
			// process symbols on next iteration, after namespace is created
			
			return true;
		};
		for(const pugi::xml_node &eDerived : dynel.by_name.getOrNew("DerivedVariable") ){
			if( !ParseDerivedVariable(log, eDerived) ) return false;
		}
		for(const pugi::xml_node &eDerived : kids.by_name.getOrNew("DerivedParameter") ){
			if( !ParseDerivedVariable(log, eDerived) ) return false; // for now
		}
		for(const pugi::xml_node &eDerived : dynel.by_name.getOrNew("ConditionalDerivedVariable") ){
			if( !ParseDerivedVariable(log, eDerived) ) return false;
		}
		
		// resolve derived variables, based on states and derived variables
		// keep symbols needed from derived variable, to create a dependency graph
		// now resolve dependencies between derived variables, to get a proper order of compoutation
		// it is really convenient, that in components without children, all namespace values other than dependent variables are readily available
		// so just look out for dependencies between derived variables of the same component type
		
		auto ResolveSymbolTable = [](const ImportLogger &log, const pugi::xml_node &eExpr, const auto &name_space, ComponentType::ResolvedTermTable &expr){
			expr.resolved.resize( expr.tab.symbol_refs.size() );
			for( int sym_seq = 0; sym_seq < (int)expr.tab.symbol_refs.size(); sym_seq++ ){
				auto sName = expr.tab.symbol_refs.at(sym_seq).c_str();
				
				auto resolved = name_space.get_id(sName);
				expr.resolved[sym_seq] = resolved;
				if( resolved < 0 ){
					log.error(eExpr, "unknown expression term %s", sName);
					return false;
				}
			}
			return true;
		};
		
		std::vector< std::set<Int> > dervar_depends_on   ( new_type.derived_variables.contents.size() );
		std::vector< std::set<Int> > dervar_depended_upon( new_type.derived_variables.contents.size() );
		for( size_t dervar_seq = 0; dervar_seq < new_type.derived_variables.contents.size(); dervar_seq++ ){
			
			auto &dervar = new_type.derived_variables.get(dervar_seq);
			
			auto ResolveAndGetDeps = [ &dervar_depends_on, &dervar_depended_upon, &dervar_seq, &ResolveSymbolTable ] ( const ImportLogger &log, const ComponentType &new_type, const auto &eDerivedValue, ComponentType::ResolvedTermTable &derived_value ) {
				if( !ResolveSymbolTable(log, eDerivedValue, new_type.name_space, derived_value) ) return false;
				
				//also gather dependencies
				for( int sym_seq = 0; sym_seq < (int)derived_value.tab.symbol_refs.size(); sym_seq++ ){	
					auto nom = new_type.name_space.get(derived_value.resolved[sym_seq]);
					if( nom.type == ComponentType::NamespaceThing::DERIVED ){
						dervar_depends_on   [ dervar_seq ].insert( nom.seq );
						dervar_depended_upon[ nom.seq ].insert( dervar_seq );
					}
					else{
						// if not derived it is known directly, skip
					}
				}
				
				//printf("%zd\n\n", derived_value.resolved.size());
				return true;
			};
			
			if( dervar.type == ComponentType::DerivedVariable::Type::VALUE ){
				if( !ResolveAndGetDeps( log, new_type, derived_nodes.get(dervar_seq), dervar.value ) ) return false;
			}
			else if( dervar.type == ComponentType::DerivedVariable::Type::CONDITIONAL ){
				auto eCase = derived_nodes.get(dervar_seq).child("Case"); // for tracking, TODO something more mantainable
				for( size_t i = 0; i < dervar.cases.size(); i++ ){
					auto &deri_case = dervar.cases[i];
					
					if( !ResolveAndGetDeps( log, new_type, eCase, deri_case.condition ) ) return false;
					if( !ResolveAndGetDeps( log, new_type, eCase, deri_case.value ) ) return false;
					eCase = eCase.next_sibling("Case");
				}
			}
			else{
				log.error(derived_nodes.get(dervar_seq), "internal error: unknown resolve dervar type");
				return false;
			}
		}
		if(log.debug_printing){
		printf("Derived variables: %zd\n", new_type.derived_variables.contents.size());
		for( size_t dervar_seq = 0; dervar_seq < new_type.derived_variables.contents.size(); dervar_seq++ ){
			printf("\t%s:", new_type.derived_variables.getName(dervar_seq) );
			for( auto depid : dervar_depends_on[dervar_seq] ){
				printf(" %s", new_type.derived_variables.getName(depid) );
			}
			printf("\n");
		}
		}
		// printf("ok\n");
		std::vector<Int> dervar_toposort, dervar_topofree;
		for( size_t dervar_seq = 0; dervar_seq < new_type.derived_variables.contents.size(); dervar_seq++ ){
			if( dervar_depends_on[dervar_seq].empty() ){
				if(log.debug_printing) printf("Free dependency: %s\n", new_type.derived_variables.getName(dervar_seq) );
				dervar_topofree.push_back(dervar_seq);
			}
		}
		
		while( !dervar_topofree.empty() ){
			auto pickout = *dervar_topofree.rbegin();
			dervar_topofree.pop_back();
			dervar_toposort.push_back(pickout);
			for( auto undep : dervar_depended_upon[pickout] ){
				dervar_depends_on[undep].erase(pickout);
				if( dervar_depends_on[undep].empty() ){
					dervar_topofree.push_back(undep);
				}
			}
		}
		
		for( size_t dervar_seq = 0; dervar_seq < new_type.derived_variables.contents.size(); dervar_seq++ ){
			if( !dervar_depends_on[dervar_seq].empty() ){
				std::string cycle_str = new_type.derived_variables.getName(dervar_seq);
				auto next_seq = dervar_seq;
				do{
					next_seq = *( dervar_depends_on[dervar_seq].begin() );
					cycle_str += " -> "; cycle_str +=  new_type.derived_variables.getName(next_seq);
				} while( next_seq != dervar_seq);
				
				log.error(derived_nodes.get(dervar_seq), "Circular dependency: %s", cycle_str.c_str());
				return false;
			}
		}
		assert( dervar_toposort.size() == new_type.derived_variables.contents.size() );
		new_type.derived_variables_topological_order = dervar_toposort;
		
		// resolve dynamics and conditions, based on states and resolved derived variables
		auto ParseStateAssignment = [&ResolveSymbolTable](const ImportLogger &log, const pugi::xml_node &eAssign, const auto &new_type, ComponentType::StateAssignment &new_assignment){
			
			auto variable = RequiredAttribute(log, eAssign, "variable");
			if(!variable) return false;
			auto value = RequiredAttribute(log, eAssign, "value");
			if(!value) return false;
			
			auto state_seq = new_type.state_variables.get_id(variable);
			if(state_seq < 0){
				log.error(eAssign, "unknown state variable %s", variable);
				return false;
			}
			// auto &state_variable = new_type.state_variables.get(state_seq);
			
			new_assignment.state_seq = state_seq;
			
			if( !ParseLemsExpression( value, new_assignment.value.tab ) ){
				log.error(eAssign, "could not parse value expression");
				return false;
			}
			// TODO non-boolean checking !
			// TODO dimension checking !
			
			// resolve right away, since namespace is known
			if( !ResolveSymbolTable(log, eAssign, new_type.name_space, new_assignment.value) ) return false;
			
			return true;
		};
		
		auto ParseDoStuff = [ &ParseStateAssignment ](const ImportLogger &log, const pugi::xml_node &eDoit, const auto &new_type, auto &do_stuff){
			
			for( const pugi::xml_node &eAssign: eDoit.children("StateAssignment") ){
				ComponentType::StateAssignment new_assignment;
				if( !ParseStateAssignment(log, eAssign, new_type, new_assignment) ) return false;
				do_stuff.assign.push_back(new_assignment);
			}
			for( const pugi::xml_node &eEvout: eDoit.children("EventOut") ){
				// send spikes
				ComponentType::EventOut out;
				out.port_seq = -1;
				
				auto port = RequiredAttribute( log, eEvout, "port" );
				if( !port ) return false;
				
				out.port_seq = new_type.event_outputs.get_id(port);
				if( out.port_seq < 0 ){
					log.error(eEvout, "unknown output event port \"%s\"", port);
					return false;
				}
				
				do_stuff.event_out.push_back(out);
			}
			for( const pugi::xml_node &eProp: eDoit.children("Transition") ){
				// LATER when regimes are to be used
				log.error(eProp, "regime Transition not supported yet");
				return false;
			}
			
			return true; // yay!
		};
		
		for(const pugi::xml_node &eDerivative : dynel.by_name.getOrNew("TimeDerivative") ){
			
			auto new_derivative = new_type.derived_variables.NewContent();
			
			auto variable = RequiredAttribute(log, eDerivative, "variable");
			if(!variable) return false;
			auto value = RequiredAttribute(log, eDerivative, "value");
			if(!value) return false;
			
			auto state_seq = new_type.state_variables.get_id(variable);
			if(state_seq < 0){
				log.error(eDerivative, "unknown state variable %s", variable);
				return false;
			}
			auto &state_variable = new_type.state_variables.get(state_seq);
			
			if(state_variable.dynamics != ComponentType::StateVariable::DYNAMICS_NONE){
				log.error(eDerivative, "time derivative for state variable %s already defined", variable);
				return false;
			}
			
			if( !ParseLemsExpression( value, state_variable.derivative.tab ) ){
				log.error(eDerivative, "could not parse value expression");
				return false;
			}
			// TODO non-boolean checking !
			// TODO dimension checking !
			
			// resolve right away, since namespace is known
			if( !ResolveSymbolTable(log, eDerivative, new_type.name_space, state_variable.derivative) ) return false;
			
			state_variable.dynamics = ComponentType::StateVariable::DYNAMICS_CONTINUOUS;
			
		}
		for(const pugi::xml_node &eOnev : dynel.by_name.getOrNew("OnEvent") ){
			
			ComponentType::OnEvent new_onevent;
			
			auto port = RequiredAttribute( log, eOnev, "port" );
			if( !port ) return false;
			
			new_onevent.in_port_seq = new_type.event_inputs.get_id(port);
			if( new_onevent.in_port_seq < 0 ){
				log.error(eOnev, "unknown input event port \"%s\"", port);
				return false;
			}
			
			if( !ParseDoStuff(log, eOnev, new_type, new_onevent) ) return false;
			
			new_type.on_events.push_back(new_onevent);
		}
		for(const pugi::xml_node &eOnco : dynel.by_name.getOrNew("OnCondition") ){
			
			ComponentType::OnCondition new_oncondition;
			
			auto test = RequiredAttribute(log, eOnco, "test");
			if(!test) return false;
			
			if( !ParseLemsExpression( test, new_oncondition.test.tab ) ){
				log.error(eOnco, "could not parse test expression");
				return false;
			}
			
			if( !new_oncondition.test.tab.terms.at(new_oncondition.test.tab.expression_root).isBoolean() ){
				log.error(eOnco, "test should be a boolean expression");
				return false;
			}
			
			// resolve right away, since namespace is known
			if( !ResolveSymbolTable(log, eOnco, new_type.name_space, new_oncondition.test) ) return false;
			
			if( !ParseDoStuff(log, eOnco, new_type, new_oncondition) ) return false;
			
			new_type.on_conditions.push_back(new_oncondition);
		}
		for(const pugi::xml_node &eInit : dynel.by_name.getOrNew("OnStart") ){
			
			for( const pugi::xml_node &eAssign: eInit.children("StateAssignment") ){
				ComponentType::StateAssignment new_assignment;
				if( !ParseStateAssignment(log, eAssign, new_type, new_assignment) ) return false;
				new_type.on_start.push_back(new_assignment);
			}
			
		}
		
		for(const pugi::xml_node &eProp : dynel.by_name.getOrNew("Regime") ){
			// TODO not used in NeuroML DB, used in core IafRef cells
			log.error(eProp, "Regime not supported yet");
			return false;
		}
		for(const pugi::xml_node &eProp : dynel.by_name.getOrNew("KineticScheme") ){
			// LATER not used in NeuroML DB, gateKS is used instead
			log.error(eProp, "KineticScheme not supported yet");
			return false;
		}
		
		// exposure is actually handled by attribute tags, just check here in case a stated Exposure is not resolved? be more flexible, to allow partial-interface odd cases being allowed whenever possible
		for(const pugi::xml_node &eProp : kids.by_name.getOrNew("Exposure") ){
			// TODO it's not actually needed anywhere, perhaps check for consistency?? LATER
			(void) eProp;
		}
		
		// check for possible exposure native names
		for(const auto &keyval : new_type.exposures.names){
			
			auto exptype_it = ComponentType::CommonExposures::names.find(keyval.first);
			if(exptype_it != ComponentType::CommonExposures::names.end()){
				auto editme = exptype_it->second;
				
				new_type.common_exposures.*editme = keyval.second;
			}
			
		}
		// check for possible requirement native names, too! (such as iCa which is not defined in concentrationModel)
		// add some type-based requirements here to avoid resolving them in the backend engine
		for(const auto &keyval : new_type.requirements.names){
			
			auto exptype_it = ComponentType::CommonRequirements::names.find(keyval.first);
			if(exptype_it != ComponentType::CommonRequirements::names.end()){
				auto editme = exptype_it->second;
				
				new_type.common_requirements.*editme = keyval.second;
			}
			
		}
		// check for common(ie natively used) names in each part of the interface
		auto GetCommonInterface = []( const auto &iface_domain, auto &iface_common_refs ){
			for(const auto &keyval : iface_domain.names){
				// a reference type is not easy to get rid of, it seems
				const auto &name_map = std::remove_reference< decltype(iface_common_refs) >::type::names;
				auto exptype_it = name_map.find(keyval.first);
				if( exptype_it != name_map.end() ){
					auto editme = exptype_it->second;
					
					iface_common_refs.*editme = keyval.second;
				}
				
			}
		};
		GetCommonInterface( new_type.event_inputs , new_type.common_event_inputs  );
		GetCommonInterface( new_type.event_outputs, new_type.common_event_outputs );
		// same for possible event inputs
		// for(const auto &keyval : new_type.event_inputs.names){
		// 	auto exptype_it = ComponentType::CommonEventInputs::names.find(keyval.first);
		// 	if(exptype_it != ComponentType::CommonEventInputs::names.end()){
		// 		auto editme = exptype_it->second;
		// 		new_type.common_event_inputs.*editme = keyval.second;
		// 	}
		// }
		
		if( log.debug_printing ){
		printf("Component type %s: \n", name);
		new_type.debug_print(dimensions);
		}
		//add type to component types! yay!
		component_types.add(new_type, name);
		
		if( log.debug_printing ) printf("Done parsing component type %s\n", name);
		
		return true;
	}
	
	bool ParseStandaloneComponentInstance( const ImportLogger &log, const pugi::xml_node &eInstance, const char *type ){
		
		Int comp_id = ParseComponentInstanceType(log, eInstance, component_types, type);
		if( comp_id < 0 ) return false;
		// perhaps run some checking on this prototype, to avoid having to run checking on each and every instance? LATER
		//printf("Tagified %s\n",type);
		// check component type, because the component instantiation information in these tags coexists with the core NeuroML properties, yay!
		//printf("Tagified %zd\n",comp_id);
		const auto &comp = component_types.get(comp_id);
		auto extends = comp.extends;
		//printf("Tagified %zd %d\n",comp_id, extends);
		if( extends == ComponentType::PURE ){
			
			// parse a general-use LEMS component prototype
			ComponentInstance new_instance;
			auto sProtoName = RequiredNmlId(log, eInstance);
			if(!sProtoName) return false;
			if( !ParseComponentInstance(log, eInstance, component_types, dimensions, type, new_instance) ) return false;
		
			// add its generic form
			if(component_instances.has(sProtoName)){
				log.error(eInstance, "standalone instance %s already defined", sProtoName);
					return false;
			}
			
			// just leave it to standalone instance prototypes, to be fully validated at instantiation time
			component_instances.add(new_instance, sProtoName);
			
			return true;
		}
		else if( extends == ComponentType::CONCENTRATION_MODEL ){
			
			return ParseConcentrationModel(log, eInstance);
		}
		else if( extends == ComponentType::ION_CHANNEL ){
			
			return ParseIonChannel(log, eInstance);
		}
		else if( extends == ComponentType::INPUT ){
			return ParseInputSource(log, eInstance);
		}
		else if( extends == ComponentType::SYNAPTIC_COMPONENT ){
			
			return ParseSynapticComponent(log, eInstance);
		}
		else if( extends == ComponentType::CELL ){
			
			return ParseArtificialCell(log, eInstance);
		}
		else{
			// no gate rates      (they are inline)
			// no gates           (they are inline)
			// no scaling factors (they are inline)
			log.error(eInstance, "standalone instances not supported for this component type");
			return false;
		}
		
	}
	
	// the constructor: attach import state to Model being imported
	ImportState(Model &_model)
		: model               ( _model                     )
		, dimensions          ( _model.dimensions          )
		, component_types     ( _model.component_types     )
		, component_instances ( _model.component_instances )
		, morphologies        ( _model.morphologies        )
		, biophysics          ( _model.biophysics          )
		, ion_species         ( _model.ion_species         )
		, conc_models         ( _model.conc_models         )
		, ion_channels        ( _model.ion_channels        )
		, cell_types          ( _model.cell_types          )
		, synaptic_components ( _model.synaptic_components )
		, input_sources       ( _model.input_sources       )
		, networks            ( _model.networks            )
		, simulations         ( _model.simulations         )
		, target_simulation   ( _model.target_simulation   )
	{
		
	}
};

//Top-level NeuroML import routine
// TODO add debug mode
bool ReadNeuroML(const char *top_level_filename, Model &model, bool entire_simulation, NmlImportContext &import_context, bool verbose, FILE *info_log, FILE *error_log){
	
	bool ok = false;
	fprintf(info_log, "Starting import from NeuroML file %s\n", top_level_filename);
	
	/* On required parse ordering:
	Cell is defined by two properties: Morphology and Biophysical. Both are required and unique in the Cell.
	Morphology explicitly specifies the 3D structure of the neuron. It is self-contained, requiring no references to other model parts.
	Biophysical property definitions may refer to specific points on the Morphology on which they are applicable g. soma only).
	Thus Morphology has to be known before interpreting that cell's Biophysical contents.
	At the same time, Morphologies, Biophysics, Channels etc. can be specified outside a particular Cell.
	The suggested order elements should appear in a NeuroML file ensures
	single files can be fully interpreted in a single-pass parse,
	although the ordering could stop holding when includes are put in the mix.
	
	Named(aka Standalone) entities shold be known across all files.
		(theoretically each file should know its include subtree, but this is both easier to code and laxer->more convenient)
	*/
	
	// all sorts of core NeuroML inputs
	static const char * known_input_types[] = {
		"pulseGenerator",
		"pulseGeneratorDL",
		"sineGenerator",
		"sineGeneratorDL",
		"rampGenerator",
		"rampGeneratorDL",
		
		// "compoundInput", TODO
		"voltageClamp",
		"voltageClampTriple",

		"timedSynapticInput",
		"poissonFiringSynapse",
		"transientPoissonFiringSynapse",
		
		"spikeArray",		
		"spikeGenerator",
		"spikeGeneratorRandom",
		"spikeGeneratorPoisson",
		"spikeGeneratorRefPoisson",
		
		"SpikeSourcePoisson",
	};
	
	static const char * known_synapse_types[] = {
		"alphaCurrentSynapse",
		"alphaSynapse",
		"expOneSynapse",
		"expTwoSynapse",
		"expThreeSynapse",
		// "doubleSynapse", TODO
		"blockingPlasticSynapse",
		
		"gapJunction",
		"silentSynapse",
		"linearGradedSynapse",
		"gradedSynapse",
		
		"expCondSynapse",
		"alphaCondSynapse",
		"expCurrSynapse",
		"alphaCurrSynapse",
	};
	
	static const char * known_artificial_cell_types[] = {
		"iafCell",
		"iafTauCell",
		"iafRefCell",
		"izhikevichCell",
		"izhikevich2007Cell",
		"adExIaFCell",
		"fitzHughNagumo1969Cell",
		"fitzHughNagumoCell",
		"pinskyRinzelCA3Cell",
		
		"IF_curr_alpha",
		"IF_curr_exp",
		"IF_cond_alpha",
		"IF_cond_exp",
		"EIF_cond_exp_isfa_ista",
		"EIF_cond_alpha_isfa_ista",
		"HH_cond_exp",
	};
	
	static const char * misc_known_top_level_types[] = {
		"notes", // not used
		"annotation", // not used
		// "property", // not used, not seen either
		"include",
		"Include",
		// "extracellularProperties", ?? but they don't seem to be referred to as id attributes?
		// "intracellularProperties", ?? but they don't seem to be referred to as id attributes?
		"morphology",
		// "IonChannelScalable", // probably not used without ion species, conductance
		"ionChannel",
		"ionChannelHH",
		"ionChannelKS",
		"concentrationModel",
		"decayingPoolConcentrationModel",
		"fixedFactorConcentrationModel",
		
		"biophysicalProperties",
		"biophysicalProperties2CaPools",
		"cell",
		"cell2CaPools",
		
		"network",
		"networkWithTemperature",
		"ComponentType",
		"Simulation",
		"Component", // will be sorted out during known & unknown component instance resolution
		"Target"
	};
	
	std::vector< const char * > known_top_level_types;
	
	AddToVec( known_top_level_types, known_input_types );
	AddToVec( known_top_level_types, known_synapse_types );
	AddToVec( known_top_level_types, known_artificial_cell_types );
	AddToVec( known_top_level_types, misc_known_top_level_types );
	
	//auxiliaries for top-level parsing
	typedef std::vector<pugi::xml_node> NodeList;
	CollectionWithNames<NodeList> top_level_nodes_by_name;
	for(auto type_name: known_top_level_types){
		top_level_nodes_by_name.add( NodeList(), type_name); //so they are always present
	}
	std::set<std::string> unknown_types;
	
	//The set of useful information retrieved
	ImportState import_state(model);
	
	// The set of files currently opened
	// import_context.clear();
	ImportLogger log(import_context, error_log, verbose); // TODO customise
	
	// The set of files to be read
	// ( since things may refer to other things referred to further along in the same or another file,
	// there is no problem with mixing the order up )
	const int MAX_INCLUDE_RECURSION = 99; // ask the devs for more, whenever needed
	
	struct FileToRead{
		std::string path;
		int include_level;
	};
	std::vector<FileToRead> files_to_read = {
		{top_level_filename,0},
		{LEMS_CoreComponents_filename,0},// see other options LATER
	};
	std::set<std::string> files_considered; // to avoid opening same path twice
	while(!files_to_read.empty()){
		
		FileToRead file_to_read = *files_to_read.rbegin();
		files_to_read.pop_back();
		const char *filename = file_to_read.path.c_str(); //obviously
		
		// TODO hide internal files ?
		fprintf(info_log, "Loading file %s\n", filename);
		//make an entry
		import_context.documents_opened.push_back( NmlFileContext() );
		NmlFileContext &cx = *import_context.documents_opened.rbegin();
		cx.filename = file_to_read.path;
		
		//load the file through Document Class
		//	NB: File data must remain resident until while the DOM exists,
		//	since all Strings returned by DOM getters are actually located in the file data buffer.
		pugi::xml_document &doc = cx.doc;
		pugi::xml_parse_result result;
		if( strcmp( filename, LEMS_CoreComponents_filename ) == 0 ){
			// the magic included file
			
			// TODO make the name directory-independent with LD script https://stackoverflow.com/a/44288493
			// extern char _binary_eden_neuroml_LEMS_CoreComponents_inc_xml_end, _binary_eden_neuroml_LEMS_CoreComponents_inc_xml_start;
			// const char * LEMS_CoreComponents_buf = &_binary_eden_neuroml_LEMS_CoreComponents_inc_xml_start;
			// size_t LEMS_CoreComponents_size = &_binary_eden_neuroml_LEMS_CoreComponents_inc_xml_end - &_binary_eden_neuroml_LEMS_CoreComponents_inc_xml_start;
			
			extern unsigned char eden_neuroml_LEMS_CoreComponents_inc_xml[]; const unsigned char *LEMS_CoreComponents_buf = eden_neuroml_LEMS_CoreComponents_inc_xml;
			extern unsigned int eden_neuroml_LEMS_CoreComponents_inc_xml_len; size_t LEMS_CoreComponents_size = eden_neuroml_LEMS_CoreComponents_inc_xml_len;
			// printf("buf %p size %zd\nload it\n", LEMS_CoreComponents_buf, LEMS_CoreComponents_size);
			// fwrite( LEMS_CoreComponents_buf, LEMS_CoreComponents_size, 1, stdout );
			result =  doc.load_buffer(LEMS_CoreComponents_buf, LEMS_CoreComponents_size); //auto management of file data
		}
		else{
			result =  doc.load_file(filename); //auto management of file data
		}
		if(!result){
			ReportErrorInFile(error_log, filename, result.offset, "Could not parse root NeuroML file: %s\n", result.description() );
			
			goto CLEANUP;
		}
		
		pugi::xml_node eRoot = doc.child("neuroml");
		if(!eRoot){
			eRoot = doc.child("Lems");
			if(!eRoot){
				fprintf(error_log, "file %s: no <neuroml> or <Lems> root element found\n", filename);
				goto CLEANUP;
			}
		}
		if(eRoot.next_sibling("neuroml") || eRoot.next_sibling("Lems")){
			fprintf(error_log, "file %s: "
				"\xF0\x9D\x96\x8D\xF0\x9D\x96\x94\xF0\x9D\x96\x89\xF0\x9D\x96\x8E\xF0\x9D\x96\x8A "
				"\xF0\x9D\x96\x93\xF0\x9D\x96\x86\xF0\x9D\x96\x99\xF0\x9D\x96\x9A\xF0\x9D\x96\x98 "
				"\xF0\x9D\x96\x8A\xF0\x9D\x96\x98\xF0\x9D\x96\x99 "
				"\xF0\x9D\x96\x97\xF0\x9D\x96\x86\xF0\x9D\x96\x89\xF0\x9D\x96\x8E\xF0\x9D\x96\x88\xF0\x9D\x96\x8E "
				"\xF0\x9D\x96\x8B\xF0\x9D\x96\x97\xF0\x9D\x96\x86\xF0\x9D\x96\x99\xF0\x9D\x96\x8A\xF0\x9D\x96\x97\n"
			, filename);
			goto CLEANUP;
		}
		
		//add backlink to entry
		void * key = eRoot.root().internal_object();
		//printf("Adding file %s with key %p\n",filename,key);
		import_context.files_by_root_element.insert(std::make_pair(key, &cx));
		
		//scan top-level elements
		for (const pugi::xml_node &eTop: eRoot.children()){
			const char *node_type = eTop.name();
			
			// check for special tag <Component> which may hide any other type
			if(strcmp(node_type, "Component") == 0){
				const char *sType = eTop.attribute("type").value();
				// is it a core NeuroML top level element ?
				// who knows, handle it anyway
				if( *sType ){
					// pass it on, NB handler should have "type" detection if they use "name"
					node_type = sType;
				}
				// else, go on with "Component"
			}

			if( !top_level_nodes_by_name.has(node_type) ){
				if( !*node_type ){
					log.error(eTop, "xml tag name missing");
					goto CLEANUP;
				}
				unknown_types.insert( std::string(node_type) );
			}
			
			// check for special tag <Include>
			if(stricmp(node_type, "include") == 0){
				const char *href = eTop.attribute("href").value();
				if(!*href){
					href = eTop.attribute("file").value();
					if(!*href){
						if( eTop.attribute("url").value() ){
								log.error(eTop,"<%s> from URL not supported yet", node_type);
						}
						else{
							log.error(eTop,"<%s> tag lacks href, file, or url attribute", node_type);
						}
						
						goto CLEANUP;
					}
				}
				
				//recursion check
				if( file_to_read.include_level >= MAX_INCLUDE_RECURSION){
					log.error(eTop,"too deep <include> recursion (reached %d)", MAX_INCLUDE_RECURSION);
					goto CLEANUP;
				}
				
				//core files check
				static const NameIndexer core_filenames = {
					{"Cells.xml", 1},
					{"Channels.xml", 2},
					{"Inputs.xml", 3},
					{"Networks.xml", 4},
					{"NeuroML2CoreTypes.xml", 5},
					{"NeuroML2CoreCompTypes.xml", 6},
					{"NeuroML2CoreDimensions.xml", 7},
					{"PyNN.xml", 8},
					{"Simulation.xml", 9},
					{"Synapses.xml", 10},
					{"NeuroMLCoreDimensions.xml", 11}, // yeah I know
				};
				
				// href needs some preprocessing before matching to core filenames, because people are funny
				const char *true_ref = href;
				// don't skip leading slashes, absolute paths must be absolute
				// skip repeated (./)* scanning from start (thus avoiding ../ or other cases)
				while( *true_ref == '.' && *true_ref == '/' ) true_ref += 2;
				// don't strip leading backslashes, because they're important on Windows (UNC paths, SMB mounts etc.)
				// fortunately people haven't spuriously added them yet
				
				if(core_filenames.count(true_ref) > 0){
					fprintf(info_log, "Skipping file %s\n", true_ref);
					// pass 
				}
				else{
					// there are two cases here, either the path is absolute (in which case it is used)
					// or the path is relative (in which case it is appended to dirname of the file that contained the <Include> tag)
					std::string resolved_path = GetRelativeFilePath(file_to_read.path, true_ref);
					
					// LATER keep path, but print non-resolved one, to avoid leaking irrelevant and possibly sensitive absolute path?
					
					// just an early check that file exists, so that the file that included the missing file can be blamed
					// don't trust the file to keep on (not) existing when it is actually opened, of course
					struct stat stat_dummy;
					if( stat( resolved_path.c_str(), &stat_dummy) < 0 ){
						log.error(eTop,"error opening file %s : %s", resolved_path.c_str(), strerror(errno) );
						goto CLEANUP;
					}
					// TODO use Windows API, instead of relying on the generosity of MinGW to contain sys/stat.h
					
					if(files_considered.count(resolved_path) > 0){
						// reopen check
						// pass
					}
					else{
						// and add new file to read
						FileToRead new_file_to_read = {resolved_path, file_to_read.include_level + 1};
						files_to_read.emplace_back(new_file_to_read);
						files_considered.insert(resolved_path);
					}
					
				}
			}
			
			top_level_nodes_by_name.getOrNew(node_type).push_back(eTop);
		}
	}
	
	// now decipher the XML elements, class by class
	// to achieve (mostly) topological order, yay!!
	
	// LEMS component types
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("ComponentType") ){
		// temp trick to skip logging for internal core components
		// and also force, TODO something more elegant...
		bool old_debug = log.debug_printing;
		bool is_internal_element = ( strcmp( log.GetFilenameFromElement(eElm), LEMS_CoreComponents_filename ) == 0 );
		if( is_internal_element ) log.debug_printing = false;
		//TODO be more selective on what to be verbose on, eg components only...
		
		
		if( !import_state.ParseLemsComponentType(log, eElm) ) goto CLEANUP;
		
		if( is_internal_element ) log.debug_printing = old_debug;
		
	}
	
	printf("Parsed all component types\n");
	printf("Parsing standalone LEMS components...\n");
	// now take one more look at unknown tags, because they might be LEMS components masquerading as standalone tags :D
	{
		
	std::set<std::string> tagified_components;
	for(const std::string &type : unknown_types){
		auto sType = type.c_str();
		// printf("Tagified %s\n",sType);
		auto comptype_id = model.component_types.get_id(sType);
		// <Component> tag included here
		// Note that special-use types like <Simulation> and <Network> are handled separately
		if( comptype_id >= 0 || type == "Component" ){
			tagified_components.insert(type);

			for( auto eElm : top_level_nodes_by_name.get(sType) ){
				if( !import_state.ParseStandaloneComponentInstance(log, eElm, sType ) ) goto CLEANUP;
			}

			
		}
		else{
			// ignore
		}
		//printf("Done tagified %s\n",sType);
	}
	// and also remove them after scanning
	for(const std::string &name : tagified_components){
		unknown_types.erase(name);
	}
	
	}
	
	// morphologies
	printf("Parsing standalone morphologies...\n");
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("morphology") ){
		if( !import_state.ParseStandaloneMorphology(log, eElm) ) goto CLEANUP; }
	
	// substance concentration models
	printf("Parsing concentration models...\n");
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("concentrationModel") ){
		if( !import_state.ParseConcentrationModel(log, eElm) ) goto CLEANUP; }
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("decayingPoolConcentrationModel") ){
		if( !import_state.ParseConcentrationModel(log, eElm) ) goto CLEANUP; }
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("fixedFactorConcentrationModel") ){
		if( !import_state.ParseConcentrationModel(log, eElm) ) goto CLEANUP; }
	
	// ion channel types
	printf("Parsing ion channels...\n");
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("ionChannel") ){
		if( !import_state.ParseIonChannel(log, eElm) ) goto CLEANUP; }
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("ionChannelHH") ){
		if( !import_state.ParseIonChannel(log, eElm) ) goto CLEANUP; }
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("ionChannelKS") ){
		if( !import_state.ParseIonChannel(log, eElm) ) goto CLEANUP; }
	
	// leave parsing for specific cell type instantiation, for biophysics
	//printf("Parsing standalone biophysics...\n");
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("biophysicalProperties") ){
		if( !import_state.ParseStandaloneBiophysics(log, eElm) ) goto CLEANUP; }
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("biophysicalProperties2CaPools") ){
		if( !import_state.ParseStandaloneBiophysics(log, eElm) ) goto CLEANUP; }
	
	// entire cell types
	printf("Parsing cell types...\n");
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("cell") ){
		if( !import_state.ParsePhysicalCell(log, eElm) ) goto CLEANUP; }
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("cell2CaPools") ){
		if( !import_state.ParsePhysicalCell(log, eElm, true) ) goto CLEANUP; }
	
	{
	// all sorts of core NeuroML artificial cells
	for(auto artificial_cell_type_name: known_artificial_cell_types){
		for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew(artificial_cell_type_name) ){
			if( !import_state.ParseArtificialCell(log, eElm) ) goto CLEANUP; }
	}
	}
	
	{
	// all sorts of core NeuroML synapses
	for(auto synapse_name: known_synapse_types){
		for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew(synapse_name) ){
			if( !import_state.ParseSynapticComponent(log, eElm) ) goto CLEANUP; }
	}
	}
	
	{
	// all sorts of core NeuroML input sources	
	for(auto input_name: known_input_types){
		for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew(input_name) ){
			if( !import_state.ParseInputSource(log, eElm) ) goto CLEANUP; }
	}
	}
	
	// networks
	printf("Parsing networks...\n");
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("network") ){
		if( !import_state.ParseNetwork(log, eElm) ) goto CLEANUP; }
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("networkWithTemperature") ){
		if( !import_state.ParseNetwork(log, eElm) ) goto CLEANUP; }
	
	// parse simulations
	printf("Parsing <Simulation>...\n");
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("Simulation") ){
		if( !import_state.ParseSimulation(log, eElm) ) goto CLEANUP; }
		
	// parse the hopefully unique target
	if( entire_simulation ){
	if(top_level_nodes_by_name.getOrNew("Target").empty()){
		fprintf(error_log, "no <Target> tag specified\n");
		goto CLEANUP;
	}
	}
	for(const pugi::xml_node &eElm : top_level_nodes_by_name.getOrNew("Target") ){
		if( !import_state.ParseTarget(log, eElm) ) goto CLEANUP; }
	
	printf("done parsing!\n");
	if(!unknown_types.empty()){
		printf("Unknown element types:\n");
		for(const std::string &name : unknown_types){
			printf("\t\"%s\"\n", name.c_str());
		}
		printf("\n");
	}
	
	//XML exploration end
	
	// TODO probe for unresolved refs in all lems (or other?) parts of the model!

	//parsing complete, no errors!
	ok = true;

	CLEANUP:
	//all handled automatically, actually
	
	return ok;
}
struct NmlImportContext_Holder::Impl : NmlImportContext { };
NmlImportContext_Holder::NmlImportContext_Holder(): impl(new Impl()){}
NmlImportContext_Holder &NmlImportContext_Holder::operator=( NmlImportContext_Holder &&rhs){
	impl=std::move(rhs.impl);
	return *this;
};
NmlImportContext_Holder::~NmlImportContext_Holder() = default;

bool ReadNeuroML(const char *top_level_filename, Model &model, bool entire_simulation, bool verbose, FILE *info_log, FILE *error_log){
	NmlImportContext import_context; // just make up one on the spot
	return ReadNeuroML(top_level_filename, model, entire_simulation, import_context, verbose, info_log, error_log);
	// beware: non-essential information like strings and xml tags, filenames etc are lost after this point, unless a holter is used!!
}
bool ReadNeuroML(const char *top_level_filename, Model &model, bool entire_simulation, NmlImportContext_Holder &import_context_holder, bool verbose, FILE *info_log, FILE *error_log){
	return ReadNeuroML(top_level_filename, model, entire_simulation, *(import_context_holder.impl), verbose, info_log, error_log);
}
