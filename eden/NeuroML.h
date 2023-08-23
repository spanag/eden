/* NeuroML API main header file */
#ifndef NEUROML_H
#define NEUROML_H

#include "Common.h"

#include "neuroml/LEMS_Expr.h"

#include <map>
#include <unordered_map>
#include <set>
#include <list>
#include <memory>

// set up a namespace LATER

// for some error reporting in the API
struct ILogProxy{
	virtual void error(const char *format, ...) const = 0;
	virtual void warning(const char *format, ...) const = 0;
	virtual ~ILogProxy(){}
};
struct NullLogProxy : public ILogProxy{
	void error(const char *format, ...) const {}
	void warning(const char *format, ...) const {}
};

//for NeuroML float and integer values
typedef float Real;
typedef long Int;

//Run-length compressed arrays.
//A list of [start, start + length) pairs.
struct IdListRle{
	std::vector<Int> run_starts;
	std::vector<Int> run_lengths;
	
	bool operator== (const IdListRle &rhs){ return run_starts == rhs.run_starts && run_lengths == rhs.run_lengths; }
	bool operator!= (const IdListRle &rhs){ return !(*this == rhs); }
	
	void Addd(Int add){
		//first index greater than add
		size_t ind = std::upper_bound(run_starts.begin(), run_starts.end(), add) - run_starts.begin(); 
		
		if(ind != run_starts.size()){
			// try merging with ind, which is greater than me
			if(add + 1 == run_starts[ind]){
				//extend downward
				run_starts[ind] = add;
				run_lengths[ind]++;
				return;
			}
		}
		
		if( ind != 0 ){
			// try merging with ind - 1, which i am greater or equal to
			int prev = ind - 1;
			Int prevend = run_starts[prev] + run_lengths[prev];
			if( prevend >= add ){
				if(prevend == add) run_lengths[prev]++;
				return;
			}
		}
		
		// otherwise extend the array :(
		run_starts.insert(run_starts.begin() + ind, add);
		run_lengths.insert(run_lengths.begin() + ind, 1);
		
		return;
	}
	void Addd(Int sta, Int len){
		//first index greater than add
		size_t ind = std::upper_bound(run_starts.begin(), run_starts.end(), sta) - run_starts.begin();
		
		if(ind != run_starts.size()){
			// try merging with ind, which is greater than me
			if(sta + len >= run_starts[ind]){
				//extend downward
				run_lengths[ind] = std::max(len, run_starts[ind] - sta + run_lengths[ind]);
				run_starts[ind] = sta;
				return;
			}
		}
		
		if( ind != 0 ){
			// try merging with ind - 1, which i am greater or equal to
			Int prev = ind - 1;
			Int prevend = run_starts[prev] + run_lengths[prev];
			if( prevend >= sta ){
				run_lengths[prev] = std::max(sta - run_starts[prev] + len, run_lengths[prev]);
				return;
			}
		}
		
		// otherwise extend the array :(
		run_starts.insert(run_starts.begin() + ind, sta);
		run_lengths.insert(run_lengths.begin() + ind, len);
		
		return;
	}
	void Add(const IdListRle &add){
		//just add the ranges one by one
		for(size_t i = 0; i < add.run_starts.size(); i++){
			Addd(add.run_starts[i], add.run_lengths[i]);
		}
		
	}
	IdListRle Minus(const IdListRle &_sub) const {
		if(_sub.run_starts.empty()) return *this; //TODO check if it needs to go after Compact()
		
		IdListRle moi = *this; moi.Compact();
		IdListRle sub = _sub; sub.Compact();
		IdListRle ret;
		
		size_t ib = 0;
		// Int sub_cover;
		for(size_t ia = 0; ia < moi.run_starts.size(); ia++){
			Int start = moi.run_starts[ia];
			Int end = moi.run_lengths[ia] + start;
			//printf("ia %zd\n", ia);
			for( ; ib < sub.run_starts.size() && sub.run_starts[ib] < end; ib++){
				//see how this interval may clip current interval
				Int cs = sub.run_starts[ib];
				Int ce = sub.run_lengths[ib] + cs;
				//printf("\tib %zd, clip %ld %ld\n", ib, start, end);
				//either: not clipping, clipping current start, splitting this interval, or fully covering it 
				if(ce < start) continue; //not clipping, proceed to next
				if(cs > start){
					//split
					//printf("\tSplit %ld, %ld - %ld\n", start, cs - start, ce);
					ret.run_starts.push_back(start);
					ret.run_lengths.push_back(cs - start);
					start = ce;
				}
				else{
					//no need to split
					//printf("\tNosplit %ld\n", ce);
				}
				start = ce; //clip start, possibly even exceeding the end (fully covered case)
				if(start >= end) break; //must exit early then, keep tis ib for next ia too
			}
			
			//also append the remaining part of the interval, if there is one
			if(start < end){
				ret.run_starts.push_back(start);
				ret.run_lengths.push_back(end - start);
			}
		}
		
		return ret;
	}
	
	//Count how many ID's are contained in the list
	//NOTE: list may be denormalized, with overlapping spans !
	Int Count() const {
		Int total = 0;
		Int counted_until = 0;
		//Sweep across the intervals with increasing interval start;
		for(size_t i = 0; i < run_starts.size(); i++){
			if(run_starts[i] + run_lengths[i] <= counted_until) continue; //already counted this entire range
			//now check if part of tis interval is already counted in total
			total += std::min(run_lengths[i], run_lengths[i] - ( counted_until - run_starts[i]) );
			counted_until = run_starts[i] + run_lengths[i]; //proceed to end of the latest input interval
		}
		return total;
	}
	void Compact(){
		if(run_starts.empty()) return;
		
		Int covered_until = run_starts[0] + run_lengths[0];
		Int new_interval_count = 0; //how many intervals have been completed in the new list?
		run_starts[new_interval_count] = run_starts[0]; //start of the first resulting segment is obvious
		for(size_t i = 1; i < run_starts.size(); i++){
			if(run_starts[i] + run_lengths[i] <= covered_until) continue; //just skip this interval
			//check if the new compact list up to this point, and the new interval, are disjoint
			if(covered_until >= run_starts[i]){
				//can merge
			}
			else{
				//they are disjoint, finalize current interval and start a new one
				run_lengths[new_interval_count] = covered_until - run_starts[new_interval_count];
				new_interval_count++;
				run_starts[new_interval_count] = run_starts[i];
			}
			//proceed to end of the latest input interval
			covered_until = run_starts[i] + run_lengths[i];
		}
		//finalize the last interval of the new compacted list
		run_lengths[new_interval_count] = covered_until - run_starts[new_interval_count];
		new_interval_count++;
		
		//and that's all the intervals there are now
		run_starts.resize(new_interval_count);
		run_lengths.resize(new_interval_count);
	}
	// write some sanity tests for these
	static bool SelfTest(){
		
		//test compacting
		IdListRle compact_test = {
			{0, 2, 3, 5, 5, 8, 9, 13, 15},
			{1, 2, 3, 2, 1, 4, 2,  1,  1}
		};
		if(compact_test.Count() != 12){
			fprintf(stderr, "compact count %s\n", std::to_string(compact_test.Count()).c_str());
			return false;
		}
		compact_test.Compact();
		const IdListRle compact_expected = {
			{0, 2, 8, 13, 15},
			{1, 5, 4,  1,  1}
		};
		if(compact_test != compact_expected){
			fprintf(stderr, "compact set %s\n", compact_test.Stringify().c_str());
			return false;
		}
		
		//test set difference
		IdListRle subtrahend = {
			{1, 4, 8, 11, 18, 21, 25, 28, 30},
			{2, 2, 1,  5,  1,  2,  2,  1,  5}
		};
		
		IdListRle minuend = {
			{3, 7, 11, 13, 15, 17, 19, 22, 25},
			{2, 3,  1,  1,  1,  1,  1,  2,  7}
		};
		
		IdListRle difference = subtrahend.Minus(minuend);
		const IdListRle difference_expected = {
			{1, 5, 12, 14, 18, 21, 32},
			{2, 1,  1,  1,  1,  1,  3}
		};
		if(difference != difference_expected){
			fprintf(stderr, "set difference %s\n", difference.Stringify().c_str());
			return false;
		}
		
		return true;
	}
	
	std::string Stringify() const {
		char tmps[100];
		std::string ret;
		for(size_t i = 0; i < run_starts.size() ; i++){
			if( run_lengths[i] == 1 ) sprintf(tmps, " %ld", run_starts[i]);
			else sprintf(tmps, " %ld-%ld", run_starts[i], run_starts[i] + run_lengths[i] - 1);
			ret += tmps;
		}
		return ret;
	}
	void debug_print(FILE *fout = stdout) const {
		fprintf(fout, "%s", Stringify().c_str());
	}
	
	// Run only in compact lists !
	bool has(Int member){
		// first index greater than member
		int ind = std::upper_bound(run_starts.begin(), run_starts.end(), member) - run_starts.begin();
		if(ind == 0){
			// all present are greater than aadd
			return false;
		}
		
		assert(run_starts[ind-1] <= member);
		Int maxi = run_starts[ind-1] + run_lengths[ind-1];
		if(member < maxi) return true;
		else return false;
	}
	template<typename Functor>
	Functor reduce(Functor doit) const {
		//Run only when list is compacted, or you may get multiple and non-monotonic indices!
		for(size_t i = 0; i < run_starts.size(); i++){
			for(Int j = 0; j < run_lengths[i]; j++){
				doit( run_starts[i] + j );
			}
		}
		return doit;
	}
	std::vector<Int> toArray() const {
		std::vector<Int> ret;
		reduce([ &ret ](Int elm){
			ret.push_back(elm);
		});
		return ret;
	}
	
	void shrink_to_fit(){
		
		run_starts.shrink_to_fit();
		run_lengths.shrink_to_fit();
	}
	
	IdListRle(){}
	IdListRle(const std::vector<Int> &_s, const std::vector<Int> _l):run_starts(_s), run_lengths(_l){}
	// make a IdListRle from a vector of integers
	template< typename T, typename = typename std::enable_if< std::is_integral<T>::value >::type > 
	IdListRle(const std::vector<T> &from_list){
		auto sorted = from_list;
		
		run_starts.assign(sorted.begin(), sorted.end());
		run_lengths.assign(sorted.size(), 1);
		std::sort(run_starts.begin(), run_starts.end());
		Compact();
		shrink_to_fit();
	}
};

//make a string-indexed, low-overhead hash table, the lazy-arcane way
struct streq{ bool operator() (const char *a, const char *b) const { return strcmp(a,b) == 0;} };
struct strhash{ size_t operator()(const char *a) const {
	//should be a string_view, pray the allocator handles this efficiently
	//should change into handmade or boost implementation, if performance warrants it
	return std::hash<std::string>()(std::string(a));
} };

template <typename Content>
using NameMap = std::unordered_map<const char *, Content, strhash, streq>;
typedef NameMap<Int> NameIndexer; //TODO replace with CollectionWithNames

template< typename Content, typename Int = long>
struct CollectionWithNames{
	typedef Content mapped_type;
	typedef NameMap<Int> _NameIndexer;
	
	std::vector<Content> contents;
	_NameIndexer names;
	std::unordered_map<Int, const char *> names_by_id;
	
	size_t size() const { return contents.size(); }
	bool empty() const { return contents.empty(); }
	bool has(Int id) const{
		return id >= 0 && id < (Int)contents.size();
	}
	bool has(const char *name) const{
		return names.count(name) > 0 ;
	}
	// this is seq really, TODO rename
	Int get_id(const char *name) const{
		if(has(name)) return names.at(name);
		else return -1;
	}
	// FIXME rename to "at" to avoid confusion with seq <-> id mappers
	const Content & get(Int id) const{
		return contents.at(id);
	}
	Content & get(Int id){
		return contents.at(id);
	}
	const Content & get(const char *name) const{
		return get(get_id(name));
	}
	Content & get(const char *name){
		return get(get_id(name));
	}
	
	// Returns empty string on none found
	const char *getName(Int id_seq) const{
		if( !has(id_seq) ) return "";
		return names_by_id.at(id_seq);
	}
	
	Int idOrNew(const char *name){
		if(has(name)) return get_id(name);
		else return add( Content(), name );
	}
	Content & getOrNew(const char *name){
		return get(idOrNew(name));
	}
	
	Int add( const Content &new_item, const char * new_name = NULL ){
		Int new_id = (Int) contents.size();
		contents.push_back(new_item);
		if(new_name){
			names.insert(std::make_pair(new_name, new_id));
			names_by_id.insert(std::make_pair(new_id, new_name));
		}
		
		return new_id;
	}
	
	// just to coax the compiler that couldn't resolve mapped_type from auto &
	Content NewContent(){return Content();}
};

// internally represented as a dense vector
// externally represented as a random set of integer ID's
template< typename Int = long>
struct BijectionToSequence{
public:
	std::vector<Int> random_ids; //ID's must be unique !
	std::unordered_map<Int, Int> sequential_by_id;
	
	Int getSequential(Int nml_id) const {
		auto it = sequential_by_id.find(nml_id);
		if(it == sequential_by_id.end()) return -1;
		return it->second;
	}
	Int getId(Int seq_id) const {
		if(!( 0 <= seq_id && seq_id < (Int)random_ids.size() )) return -1;
		else return random_ids[seq_id];
	}
	
	bool hasId(Int id) const {
		return (sequential_by_id.count(id) > 0);
	}
	
protected:
	bool add(Int new_id){
		random_ids.push_back(new_id);
		sequential_by_id.insert(std::make_pair(new_id, random_ids.size() - 1));
		return true;
	}
};

template< typename Content, typename Int = long>
struct CollectionWithIds : public BijectionToSequence<Int> {
	std::vector<Content> contents;
	
	Content &atSeq(Int seq){
		return contents[seq];
	}
	const Content &atSeq(Int seq) const {
		return contents[seq];
	}
	const Content &atId(Int id) const {
		// see http://www.open-std.org/jtc1/sc22/wg21/docs/cwg_defects.html#108
		// in short words, when the base is a template its contents may differ among explicit specializations
		// so derived cannot be sure about the base's members existing
		return  contents[ BijectionToSequence<Int>::getSequential(id)] ;
	}
	
	bool add(const Content& new_elm, Int new_id){
		contents.push_back(new_elm);
		BijectionToSequence<Int>::add(new_id);
		return true;
	}
	
	size_t size() const {
		return contents.size();
	}
	
};

//Physical quantities in use
struct Dimensionless{
	static const char *NAME;
};
struct Time{
	static const char *NAME;
};
struct Frequency{
	static const char *NAME;
};
struct Length{
	static const char *NAME;
};
struct Voltage{
	static const char *NAME;
};
struct Current{
	static const char *NAME;
};
struct Resistance{
	static const char *NAME;
};
struct Conductance{
	static const char *NAME;
};
struct Resistivity{
	static const char *NAME;
};
struct Conductivity{
	static const char *NAME;
};
struct Capacitance{
	static const char *NAME;
};
struct SpecificCapacitance{
	static const char *NAME;
};
struct Concentration{
	static const char *NAME;
};
struct Permeability{
	static const char *NAME;
};
struct RhoFactor{
	static const char *NAME;
};
struct Temperature{
	static const char *NAME;
};
// my compiler doesn't have template globals yet :(   ( it now does, fix LATER )


// Get scaling factors for various measurement units of physical quantities, to support NeuroML
struct ScaleEntry{
	const char* name; int pow_of_10; double scale; double offset;
	
	// operators work only for absolute linear scales! Do not use with non-kelvin temperatures or other weird offset units!
	ScaleEntry operator*( const ScaleEntry &rhs ) const {
		assert(offset == 0 && rhs.offset == 0);
		ScaleEntry ret = { "Derived", pow_of_10 + rhs.pow_of_10, scale * rhs.scale, 0 };
		return ret;
	}
	ScaleEntry operator/( const ScaleEntry &rhs ) const {
		assert(offset == 0 && rhs.offset == 0);
		ScaleEntry ret = { "Derived", pow_of_10 - rhs.pow_of_10, scale / rhs.scale, 0 };
		return ret;
	}
	ScaleEntry operator^( int powah ) const {
		assert(offset == 0);
		ScaleEntry ret = { "Derived", pow_of_10 * powah, std::pow(scale, powah), 0 };
		return ret;
	}
	ScaleEntry operator^( double powah ) const {
		assert(offset == 0);
		double truepow = (pow_of_10 * powah);
		int integerpow = int(truepow);
		double residue = truepow - integerpow;
		ScaleEntry ret = { "Derived", int(pow_of_10 * powah), pow10(residue) * std::pow(scale, powah), 0 };
		return ret;
	}
	// these do work for every unit, though
	ScaleEntry to( const ScaleEntry &to ) const {
		ScaleEntry ret = { "Derived", pow_of_10 - to.pow_of_10, scale / to.scale,  (offset - to.offset) / (to.scale * pow10(to.pow_of_10)) };
		return ret;
	}
	double ConvertTo( double value, const ScaleEntry &to ) const {
		return value * (scale / to.scale) * pow10(pow_of_10 - to.pow_of_10) + ( (offset - to.offset) / (to.scale * pow10(to.pow_of_10)) );
	}
}; // SI value = this value * scale * 10^pow_of_10 + offset
typedef const ScaleEntry (ScaleList)[];

template<typename UnitType> struct Scales{
	static const ScaleEntry native; //the simulator's internal unit name
	static const ScaleList scales; //supported unit names
};



//------------------> Parsed representations of NeuroML entities

// LATER forward decalrations

// NOTE: all char * strings in the Model are invalid after the XML files are parsed
// TODO: add a persistent String table for easier diagnostics

struct Base{
	/*const char *id;
	bool FillIn(const ImportLogger &log, pugi::xml_node from_node){
		id = RequiredNmlId(log, from_node);
		return !!id;
	}*/
};
struct Standalone : public Base{
	
};

struct Dimension{
	int m, l, t, i, k, n, j;
	// String[] symbols = {"kg", "m", "s", "A", "K", "mol", "cd"};
	static constexpr int Dimension::* members[] = {
		&Dimension::m,
		&Dimension::l,
		&Dimension::t,
		&Dimension::i,
		&Dimension::k,
		&Dimension::n,
		&Dimension::j
	};
	static Dimension Unity() {
		return Dimension (0, 0, 0, 0, 0, 0, 0);
	}
	std::string Stringify() const {
		std::string ret;
		
		auto AddDim = [&ret](int dim, const char *dimname){
			if(dim == 0) return;
			if(!ret.empty()){
				ret += " * ";
			}
			ret += dimname + ("^" + std::to_string(dim));
		};
		
		AddDim(m, "m");
		AddDim(l, "l");
		AddDim(t, "t");
		AddDim(i, "i");
		AddDim(k, "k");
		AddDim(n, "n");
		AddDim(j, "j");
		
		if(ret.empty()) ret = "unitless";
		return ret;
	}
	Dimension(){
		//reset to unity, perhaps invalid-value LATER?
		for( auto editme : members ){
			this->*editme = 0;
		}
	}
	Dimension( int _m, int _l, int _t, int _i, int _k, int _n, int _j ){ 
		m = _m; l = _l; t = _t; i = _i; k = _k; n = _n; j = _j;
	}
	
	bool operator== (const Dimension &rhs) const{
		for( auto editme : members ){
			if( this->*editme != rhs.*editme ) return false;
		}
		return true;
	}
	bool operator!= (const Dimension &rhs) const{
		return !( *this == rhs );
	}
	Dimension operator* (const Dimension &rhs) const{
		Dimension ret;
		for( auto editme : members ){
			ret.*editme = this->*editme + rhs.*editme;
		}
		return ret;
	}
	Dimension operator/ (const Dimension &rhs) const{
		Dimension ret;
		for( auto editme : members ){
			ret.*editme = this->*editme - rhs.*editme;
		}
		return ret;
	}
};

// TODO refactor all to parametric instead of compile types
extern const Dimension LEMS_Voltage;
extern const Dimension LEMS_Time;
extern const Dimension LEMS_Concentration;
extern const Dimension LEMS_Area;

struct LemsUnit : public ScaleEntry{
	std::string name;
	LemsUnit( std::string _name, int pow, double sca = 1.0, double off = 0.0 ){
		name = _name;
		pow_of_10 = pow;
		scale = sca;
		offset = off;
	}
	
	LemsUnit(const ScaleEntry &from)
		: ScaleEntry(from) {
		name = from.name;
		ScaleEntry::name = name.c_str(); // for good measure
	}
	LemsUnit() = delete;
};

struct DimensionSet{
	
	struct LEMS_DimensionLessThan{ bool operator() ( const Dimension &a, const Dimension &b) const {
		// lexicographic comparison
		for( auto editme : a.members ){
			if( a.*editme < b.*editme ) return true;
			if( a.*editme > b.*editme ) return false;
		}
		return false; // tie breaker
	}};
	
	struct DimensionInfo{
		// other than Dimension which is the primary key
		std::string name; // multi-name tracking for convenience LATER
		LemsUnit native;
		std::vector<LemsUnit> units;
		
		template<typename UnitType>  
		static DimensionInfo FromOld(std::string _name){
			DimensionInfo ret(_name, Scales<UnitType>::native);
			for(auto _scale : Scales<UnitType>::scales){
				ret.units.push_back(_scale);
			}
			return ret;
		}
		DimensionInfo( std::string _name )
			: name(_name),  native("SI units", 0, 1.0, 0) {
				
		}
		DimensionInfo(std::string _name, const LemsUnit &_native, const std::vector<LemsUnit> &_units = std::vector<LemsUnit>())
			: name(_name), native(_native), units(_units) {
				
		}
	};
	
	std::map< std::string, Dimension > dimensions_by_name; // more than one name may map to the same dimension
	std::map< Dimension, DimensionInfo, LEMS_DimensionLessThan > info_by_dimension;
	
	bool Add( Dimension dim, DimensionInfo info ){
		bool already_set = Has(info.name) || Has(dim);
		
		dimensions_by_name.insert(std::make_pair( info.name, dim ));
		if(!Has(dim)) info_by_dimension.insert(std::make_pair( dim, info )); // official dimension name is the first encountered
		
		return !already_set;
	}
	
	bool AddUnit(Dimension dim, LemsUnit unit){
		if(!Has(dim)) return false;
		info_by_dimension.at(dim).units.push_back(unit); // perhaps check for duplicates LATER
		return true;
	}
	
	bool Has(const std::string &name) const {
		return dimensions_by_name.count(name) > 0;
	}
	bool Has(Dimension dim) const {
		return info_by_dimension.count(dim) > 0;
	}
	
	// return dummies on missing
	// will consolidate on full info output LATER
	Dimension Get(std::string name) const {
		if(!Has(name)) return Dimension();
		return dimensions_by_name.at(name);
	}
	std::string GetName(Dimension dim) const {
		if(!Has(dim)) return "";
		return info_by_dimension.at(dim).name;
	}
	const std::vector<LemsUnit> &GetUnits( Dimension dim ) const { 
		const static auto Nothing = std::vector<LemsUnit>();
		if(!Has(dim)) return Nothing;
		return info_by_dimension.at(dim).units;
	}
	// assume fundamental units for unknown dimensions
	const LemsUnit &GetNative( Dimension dim ) const {
		if(!Has(dim)) return GetNative(Dimension::Unity());
		return info_by_dimension.at(dim).native;
	}
	
	std::string Stringify( Dimension dim ) const {
		if(Has(dim)) return GetName(dim);
		else return dim.Stringify();
	}
	
	DimensionSet(){
		AddDefaults();
	}
	
private:
	void AddDefaults();
};

struct ComponentType{
	
	enum CoreType{
		PURE = 0, // extends nothing
		
		GATE_RATE, 
		GATE_TAU, 
		GATE_INF, 
		GATE, 
		ION_CHANNEL,
		CONDUCTANCE_SCALING,
		
		CONCENTRATION_MODEL,
		
		// this is a tricky point, as both inputs, continuous and spike-based synapses all inherit from basePointCurrent
		// so detection should be through validating the component interface whenever they are used, and not a priori
		SYNAPTIC_COMPONENT,
		
		BLOCK_MECHANISM,
		PLASTICITY_MECHANISM,
		
		INPUT,
		
		CELL // an "artificial" cell, since "biophysical" ones cannot be fully captured in LEMS
	};
	
	struct EventPortIn{
		// whatever, if there is any difference between them
	};
	
	struct EventPortOut{
		// whatever, if there is any difference between them
	};
	
	struct BaseNamedProperty{
		// name is implicitly added on containers 
		Dimension dimension;
	};
	struct BaseValue{
		Real value; // to engine units, what does this mean though?
	};
	
	
	// Parameter requires the value to be specified at instantiation,
	// while Constant requires the value to be specified inline in the component definition.
	// Property can be specified both ways, instantiation value takes priority over the one specified inline.
	struct Constant : public BaseNamedProperty, public BaseValue{
		// value must be set (ie not NaN) at all times
	};
	// includes Parameter sub-case
	struct Property : public BaseNamedProperty, public BaseValue{
		// value is NaN if missing, and has to be specified at instantiation time
	};
	struct Exposure // : public BaseNamedProperty
	{
		// must map to a direct or derived variable!
		enum Type{
			NONE, // must resolve!
			STATE,
			DERIVED
			// no need to expose properties since they can be exposed through derived variables
		} type;
		Int seq; // in respective container
		
		std::string Stringify() const {
			std::string ret;
			if(type == NONE){
				ret += "none??";
			}
			else if(type == STATE){
				ret += "state";
			}
			else if(type == DERIVED){
				ret += "derived";
			}
			else{
				ret += "unknown";
			}
			ret += " ";
			ret += std::to_string(seq);
			
			return ret;
		}
	};
	struct Requirement : public BaseNamedProperty{
		// and that's all
	};
	
	struct ResolvedTermTable{
		TermTable tab;
		std::vector<Int> resolved; // TermTable symbols to namespace things
	};
	
	struct StateVariable : public BaseNamedProperty{
		
		enum Dynamics{
			DYNAMICS_NONE,
			DYNAMICS_CONTINUOUS,
		};
		
		Dynamics dynamics;
		ResolvedTermTable derivative;
		// unfortunately LEMS has no initial variable property, they will be initialized through OnStart actions
		
		StateVariable(){
			dynamics = DYNAMICS_NONE;
		}
	};
	
	struct DerivedVariable : public BaseNamedProperty{
		
		struct Case{
			ResolvedTermTable condition;
			ResolvedTermTable value;
		};
		
		enum Type{
			NONE,
			VALUE,
			CONDITIONAL,
			SELECT
		} type;
		// implement select LATER when child elements are supported
		
		ResolvedTermTable value; // if value
		std::vector<Case> cases; // if many conditions
		int default_case; // if many conditions
		// NB if no case conditions apply, then is 0 : https://github.com/LEMS/jLEMS/blob/master/src/main/java/org/lemsml/jlems/core/eval/ConditionalDBase.java#L29
		DerivedVariable() : BaseNamedProperty(){
			type = NONE;
			default_case = -1;
		}
	};
	
	struct StateAssignment{
		Int state_seq;
		ResolvedTermTable value;
	};
	struct EventOut{
		Int port_seq;
	};
	
	struct OnCondition{
		ResolvedTermTable test;
		
		std::vector<StateAssignment> assign;
		std::vector<EventOut> event_out;
	};
	
	struct OnEvent{
		Int in_port_seq;
		
		std::vector<StateAssignment> assign;
		std::vector<EventOut> event_out;
	};
	
	// the namespace of a LEMS Component consists of:
	// 	constants, parameters etc.
	// 	requirements (not exposures, these are exterior only)
	// 	state variables and derived quantities
	// 	attached Components such as children, synapses, gap junction peers etc.
	struct NamespaceThing{
		enum Type{
			NONE,
			
			// LEMS defined types of stuff
			CONSTANT,
			PROPERTY,
			REQUIREMENT,
			// EXPOSURE, do not belong to namespace
			STATE,
			DERIVED,
			
			// nonstandard EDEN extension, but it's necessary for real models, oh well
			VARREQ,
			
			// magical requirements are in CommonRequirements
			
		}type;
		static const char *getTypeName(Type type) {
			switch(type){
				case CONSTANT   : return "Constant";
				case PROPERTY   : return "Property";
				case REQUIREMENT: return "Requirement";
				case STATE      : return "StateVariable";
				case DERIVED    : return "DerivedVariable";
				default: return "Invalid";
			}
		}
		const char *getTypeName() const { return getTypeName(type); }
		Int seq; // in respective container
	};
	
	template< typename ConcreteLister >
	struct AbstractLister{
	public:
		bool Superset(const ConcreteLister &rhs ) const {
			for( auto editme : ConcreteLister::members ){
				if( ( (ConcreteLister *)this->*editme < 0 ) && !( rhs.*editme < 0 ) ) return false;
			}
			return true;
		}
		
	protected:
		AbstractLister(){
			for( auto editme : ConcreteLister::members ){
				(ConcreteLister *)this->*editme = -1;
			}
		}
	private:
		friend ConcreteLister;
	};
	
	struct CommonRequirements : public AbstractLister<CommonRequirements>{
		using AbstractLister<CommonRequirements>::AbstractLister;
		// all are requirement container seq id's
		
		// for anything
		Int time;
		Int temperature;
		
		Int membrane_voltage; // for most compartment-associated components. NB only one dimensionality can be exposed !
		Int membrane_surface_area; // for childern of <Segment>: ConcentrationModels and IonChannels
		
		Int external_current; // for artificial cells. NB only one dimensionality can be exposed !
		
		Int ion_current; // for concentration models
		Int initial_concentration_intra;
		Int initial_concentration_extra;
		
		// for many calcium-dependent, compartment-associated components
		Int calcium_concentration_intra;
		
		// for gate rates that need other gate rates to work
		Int gate_rate_alpha;
		Int gate_rate_beta ;
		
		Int gate_rate_scale; // Q10-based scale factor, for gate rates that factor it in
		
		// Int instances; // for gates, a required parameter actually
		
		Int peer_voltage; // for gap junctions
		
		Int block_factor; // for home-made plastic synapse without children
		Int plasticity_factor; // for home-made plastic synapse without children
		
		static constexpr Int CommonRequirements::* members[] = {
			&CommonRequirements::time,
			&CommonRequirements::temperature,
			&CommonRequirements::membrane_voltage,
			&CommonRequirements::membrane_surface_area,
			&CommonRequirements::external_current,
			&CommonRequirements::ion_current,
			&CommonRequirements::initial_concentration_intra,
			&CommonRequirements::initial_concentration_extra,
			&CommonRequirements::calcium_concentration_intra,		
			&CommonRequirements::gate_rate_alpha,		
			&CommonRequirements::gate_rate_beta ,		
			&CommonRequirements::gate_rate_scale,		
			&CommonRequirements::peer_voltage,		
			&CommonRequirements::block_factor,		
			&CommonRequirements::plasticity_factor,		
		};
		
		const static NameMap< Int ComponentType::CommonRequirements::* > names;
	};
	
	struct CommonExposures : public AbstractLister<CommonExposures>{
		using AbstractLister<CommonExposures>::AbstractLister;
		// all are exposure container seq id's
		
		// for artificial cells
		Int membrane_voltage;
		
		// for current-based synapstic components ( includes conductance-based models ):
		Int current;
		//Int conductance;
		
		// for whole ion channel components, and conductance-based synaptic components or inputs
		Int conductance;
		
		// for gate variable and kinetic scheme rates
		Int rate;
		Int time;
		Int variable;
		
		// for gate components
		Int gate_variable;
		// for gate and ion channel components
		Int conductivity_fraction; // of gbar
		
		// for concentration models
		Int concentration_intra;
		Int concentration_extra;
		
		// for scaling fctors (such as ion channel conductance)
		Int scaling_factor;
		
		// for block mechanisms in blocking/plastic synapses
		Int block_factor;
		
		// for block mechanisms in blocking/plastic synapses
		Int plasticity_factor;
		
		
		static constexpr Int CommonExposures::* members[] = {
			&CommonExposures::membrane_voltage,
			&CommonExposures::current, 
			&CommonExposures::rate,
			&CommonExposures::time,
			&CommonExposures::variable,
			&CommonExposures::gate_variable,
			&CommonExposures::conductivity_fraction,
			&CommonExposures::conductance,
			&CommonExposures::concentration_intra,
			&CommonExposures::concentration_extra,
			&CommonExposures::scaling_factor,
			&CommonExposures::block_factor,
			&CommonExposures::plasticity_factor,
		};
		
		const static NameMap< Int ComponentType::CommonExposures::* > names;
	};
	
	struct CommonEventInputs : public AbstractLister<CommonEventInputs>{
		using AbstractLister<CommonEventInputs>::AbstractLister;
		// all are event_in container seq id's
		
		// spikes for synapses
		Int spike_in;
		
		static constexpr Int CommonEventInputs::* members[] = {
			&CommonEventInputs::spike_in,
		};
		
		const static NameMap< Int ComponentType::CommonEventInputs::* > names;
	};
	struct CommonEventOutputs : public AbstractLister<CommonEventOutputs>{
		using AbstractLister<CommonEventOutputs>::AbstractLister;
		// all are event_out container seq id's
		
		// spikes for synapses
		Int spike_out;
		
		static constexpr Int CommonEventOutputs::* members[] = {
			&CommonEventOutputs::spike_out,
		};
		
		const static NameMap< Int ComponentType::CommonEventOutputs::* > names;
	};
	
	CoreType extends;
	
	static std::list<std::string> eternal_strings; // not necessary thanks to persistent import context, use that or make a NameTable struct LATER
	
	
	// Constants may not be defined for each instance separately.
	CollectionWithNames<Constant> constants;
	// Parameters, including <Property>s, may be defined for each component instance separately.
	CollectionWithNames<Property> properties;
	
	// Internal dynamics and "assigned variables" (in NEURON MOD file terminology)
	CollectionWithNames<StateVariable> state_variables;
	CollectionWithNames<DerivedVariable> derived_variables;
	std::vector<Int> derived_variables_topological_order; // topological order is a sequence of id's, not a mapping id -> order
	
	// nonstandard EDEN extension, but it's necessary for real models, oh well
	CollectionWithNames<Requirement> variable_requirements;
	
	CollectionWithNames<NamespaceThing> name_space;
	
	std::vector<StateAssignment> on_start;
	std::vector<OnCondition> on_conditions;
	std::vector<OnEvent> on_events;
	
	// Interface of LEMS component
	CollectionWithNames<Requirement> requirements;
	CommonRequirements common_requirements;
	CollectionWithNames<Exposure> exposures;
	CommonExposures common_exposures;
	CollectionWithNames<EventPortIn> event_inputs;
	CommonEventInputs common_event_inputs;
	CollectionWithNames<EventPortOut> event_outputs;
	CommonEventOutputs common_event_outputs;
	// same with event_outputs for artificial cells LATER
	
	
	const Dimension getExposureDimension(const char *name) const {
		if( !exposures.has(name) ) return Dimension();
		return getExposureDimension( exposures.get_id(name) );
	}
	const Dimension getExposureDimension( Int exp_seq ) const{
		// perhaps invalidtype internal error LATER?
		if( !exposures.has(exp_seq) ) return Dimension();
		auto exp = exposures.get(exp_seq);
		if(exp.type == Exposure::STATE){
			return state_variables.get(exp.seq).dimension;
		}
		else if(exp.type == Exposure::DERIVED){
			return derived_variables.get(exp.seq).dimension;
		}
		else{
			// internal error
			return Dimension();
		}
	}
	const Dimension getNamespaceEntryDimension(Int seq) const {
		// perhaps invalidtype LATER?
		if( !name_space.has(seq) ) return Dimension();
		auto exp = name_space.get(seq);
		if(exp.type == NamespaceThing::CONSTANT){
			return constants.get(exp.seq).dimension;
		}
		else if(exp.type == NamespaceThing::STATE){
			return state_variables.get(exp.seq).dimension;
		}
		else if(exp.type == NamespaceThing::DERIVED){
			return derived_variables.get(exp.seq).dimension;
		}
		else if(exp.type == NamespaceThing::PROPERTY){
			return properties.get(exp.seq).dimension;
		}
		else if(exp.type == NamespaceThing::REQUIREMENT){
			return requirements.get(exp.seq).dimension;
		}
		else if(exp.type == NamespaceThing::VARREQ){
			return variable_requirements.get(exp.seq).dimension;
		}
		else{
			assert(false);
			return Dimension();
		}
	}
	
	// dimensions of native exposures
	bool GetExposureAndDimension( Int exp_seq, Dimension &dim_out ) const;
	bool GetRequirementAndDimension( Int req_seq, Dimension &dim_out ) const {
		if( !requirements.has(req_seq) ) return false;
		else{
			dim_out = requirements.get(req_seq).dimension;
			return true;
		}
	}
	
	bool GetCurrentOutputAndDimension( Dimension &dim_out ) const;
	bool GetVoltageRequirementAndDimension( Dimension &dim_out ) const;
	
	
	void debug_print( const DimensionSet &dimensions ) const ;
};
struct ComponentInstance{
	// closure of a realized instantiation, with e.g. parameter values for the specific instance
	struct ParameterOverride{
		Int seq; // of the component's parameters
		Real value;
	};
	
	Int id_seq; // of the component
	
	std::vector<ParameterOverride> parms;
	
	bool ok() const { return ( id_seq >= 0 ); }
	void clear() { id_seq = -1; }
	
	ComponentInstance(){ clear(); }
};

// A part of a whole synapse, actually;
// The mechanism attached to the pre- or post-synaptic component
// For example, a one-way spiking synapse has only one component, 
// on the post-synaptic side; spikes are automatically sent when a spiking Connection is defined, with no additional LEMS component
// (an intenal NetCon is spawned for each such instance, see https://github.com/NeuroML/org.neuroml.export/blob/master/src/main/java/org/neuroml/export/neuron/NeuronWriter.java#L607 )
// Continuous and spiking components differ only in the ir handling of peer-> values, and perhaps event ports
struct SynapticComponent{
	
	enum Type{
		NONE = 0,
		
		SILENT,
		
		ALPHA_CURRENT,
		ALPHA,
		EXP,
		EXPTWO,
		EXPTHREE,
		
		GAP,
		GAPHALF,
		GRADED_PRINZ_ET_AL,
		BLOCKING_PLASTIC,
		
		PYNN_EXP_COND,
		PYNN_ALPHA_COND,
		PYNN_EXP_CURR,
		PYNN_ALPHA_CURR,
		
		COMPOSITE, // TODO
		COMPONENT, // whatever LEMS
		MAX // just a sentinel to count how many types exist
	} type;
	
	//used for conductance-based synapses
	struct SilentSynapse{}; //does precisely nothing, useful for one-way continuous synapses
	
	//--> conductance-based chemical synapses
	struct BaseConductanceBasedSynapse{
		Real gbase;
		Real erev;
	};
	
	// tripleexp decays with two different conductances, whatever that means
	struct BaseConductanceBasedSynapseTwo : public BaseConductanceBasedSynapse{
		Real gbase2;
	};
	
	struct AlphaSynapse : public BaseConductanceBasedSynapse{
		Real tau;
	};
	struct ExpOneSynapse : public BaseConductanceBasedSynapse{
		Real tauDecay;
	};
	struct ExpTwoSynapse : public ExpOneSynapse{
		Real tauRise;
	};
	struct ExpThreeSynapse : public ExpTwoSynapse{
		Real tauDecay2;
		Real gbase2;
	};
	struct BlockingPlasticSynapse : public ExpTwoSynapse{
		struct BlockMechanism{
			enum Type{
				NONE,
				VOLTAGE_CONC_DEP,
				COMPONENT// TODO
			} type;
			
			Int species; // for VOLTAGE_CONC_DEP
			
			Real blockConcentration; // for VOLTAGE_CONC_DEP
			Real scalingConc; // for VOLTAGE_CONC_DEP
			Real scalingVolt; // for VOLTAGE_CONC_DEP - whould there not be an offset voltage too ?
			
			ComponentInstance component; // just for lemsification, for now
			
			BlockMechanism(){ type = NONE; }
		};
		struct PlasticityMechanism{
			enum Type{
				NONE,
				TSODYKS_MARKRAM_DEP,
				TSODYKS_MARKRAM_DEP_FAC,
				COMPONENT// TODO
			} type;
			
			Real initReleaseProb; // for TSODYKS_MARKRAM_DEP and FAC
			Real tauRec; // for TSODYKS_MARKRAM_DEP and FAC
			Real tauFac; // for TSODYKS_MARKRAM_DEP_FAC
			
			ComponentInstance component; // just for lemsification, for now
			
			PlasticityMechanism(){ type = NONE; }
		};
		// TODO could there be multiple of those in a row?
		BlockMechanism block_mechanism;
		PlasticityMechanism plasticity_mechanism;
		
		BlockingPlasticSynapse(){ };
	};
	struct PynnSynapse{
		Real e_rev;
		Real tau_syn;
	};
	
	//--> graded synapses
	struct GapJunction{
		Real conductance; //without erev
	};
	
	// intended to be an one-way gap junction
	struct LinearGradedSynapse : public GapJunction{  };
	
	struct GradedSynapsePrinzEtAl : public BaseConductanceBasedSynapse{
		Real delta;
		Real k;
		Real Vth;
	};
	
	//--> pure-current synapses
	struct AlphaCurrentSynapse{
		Real ibase;
		Real tau;
	};
	
	bool HasVpeer( const CollectionWithNames<ComponentType> &component_types ) const;
	bool HasSpikeIn( const CollectionWithNames<ComponentType> &component_types ) const;
	
	bool GetCurrentOutputAndDimension( const CollectionWithNames<ComponentType> &component_types, Dimension &dim_out ) const;
	bool GetVoltageInputAndDimension( const CollectionWithNames<ComponentType> &component_types, Dimension &dim_out ) const;
	bool GetVpeerInputAndDimension( const CollectionWithNames<ComponentType> &component_types, Dimension &dim_out ) const;
	
	// if native
	union{
		AlphaSynapse alpha;
		ExpThreeSynapse exp;
		LinearGradedSynapse gap;
		GradedSynapsePrinzEtAl graded;
		AlphaCurrentSynapse alpha_current;
		PynnSynapse pynn;
	};
	// for blocking plastic
	BlockingPlasticSynapse blopla;
	
	// otherwise
	ComponentInstance component;
};

// Spiking input sources send spikes to target cell segments
// (just as if a presynaptic compartment had sent them due to voltage exceeding spikeThresh)
struct InputSource{
	
	struct Spike{
		Real time_of_occurrence;
	};
	
	enum Type{
		PULSE,// single DC pulse
		PULSE_DL, // same but dimensionless
		SINE, // a finite-duration sine burst
		SINE_DL, // same but dimensionless		
		RAMP, // steady DC current with a finite-duration ramp part
		RAMP_DL, // same but dimensionless	
		VOLTAGE_CLAMP, // clamps (non-ideally) for a finite duration
		VOLTAGE_CLAMP_TRIPLE, // clamps (non-ideally) at three different voltgae levels
		TIMED_SYNAPTIC, // a spike list with attached synapse that the spikes trigger
		POISSON_SYNAPSE, // a spike generator with attached synapse that it triggers
		POISSON_SYNAPSE_TRANSIENT, // like above but within a fixed time window
		
		SPIKE_LIST, // list of spikes in time, sort of like timedSynapticInput
		SPIKE_PERIODIC, // periodic repetition of spikes
		SPIKE_RANDOM, // spikes in uniformly random intervals
		SPIKE_POISSON, // spikes in poisson random intervals
		SPIKE_POISSON_REF, // spikes in ((poisson random) + min_delay) intervals
		
		PYNN_SPIKE_POISSON, // a finite-duration burst of spikes in poisson random intervals
		
		COMPOSITE, // TODO
		COMPOSITE_DL, // TODO
		COMPONENT,
		MAX // just a sentinel to count how many types exist
	} type;
	
	Real amplitude; // or baseline for ramp input
	
	Real duration; // for finite-duration events
	Real delay; // for finite-duration events
	Real period; // for sine and periodic spike
	
	Real phase; // for sine
	Real startAmplitude; // for ramp
	Real finishAmplitude; // for ramp
	
	Real active; // for triple voltage clamp
	Real conditioningVoltage; // for triple voltage clamp
	Real testingVoltage; // for simple and triple voltage clamp
	Real returnVoltage; // for triple voltage clamp
	Real seriesResistance; // for simple and triple voltage clamp
	
	Real averageRate; // for poisson spike inputs
	Real maxISI, minISI;// for bounded ISI inputs
	
	Int synapse; // synaptic component, for firing-synapse sorts of input. Component should accept spikes to work.
	
	std::vector<Spike> spikes; // for spike-list events. Must be kept sorted to increasing time.
	
	bool GetVoltageInputAndDimension( const CollectionWithNames<ComponentType> &component_types, const CollectionWithNames<SynapticComponent> &synapses, Dimension &dim_out ) const;
	bool GetCurrentOutputAndDimension( const CollectionWithNames<ComponentType> &component_types, const CollectionWithNames<SynapticComponent> &synapses, Dimension &dim_out ) const;
	bool HasSpikeOut( const CollectionWithNames<ComponentType> &component_types ) const;
	
	ComponentInstance component; // for LEMS inputs
};

// actually NeuroML holds no information about the ions, except for name?
struct IonSpecies{
	// This is a dummy class just in case NeuroML adds an ion library or something
};

// The shape of a physical cell
struct Morphology{
	
	struct Segment{
		struct Point3DWithDiam{
			Real x,y,z,d;
		};
		//Int id; not really necessary
		Int parent; // sequential index, not id !
		Real fractionAlong;
		Point3DWithDiam proximal, distal;
		
		void debug_print() const {
			printf("parent: %3ld proximal:(%3.3f, %3.3f, %3.3f), %2.3f distal:(%3.3f, %3.3f, %3.3f), %2.3f ",
				parent, proximal.x, proximal.y, proximal.z, proximal.d,
				distal.x, distal.y, distal.z, distal.d
			);
		}
	};
	
	// Segment ID's are increasing but may be non-sequential !
	// This means NeuroML ID's must be mapped to position in the segment array !
	// All internal references are resolved serial positions in the segment list, to preserve sanity
	//std::vector<Segment> segment_list;
	//BijectionToSequence ids_vs_segs;
	
	CollectionWithIds<Segment> segments;
	
	struct InhomogeneousParameter{
		// Note: This concept is a direct port of the SubsetDomainIterator class and its style parameters, found in  NEURON's CellBuilder scripts.
		// SubsetDomainIterator Documentation: share/lib/hoc/celbild/celset.hoc , proc hints() (NEURON source tree)
		// SubsetDomainIterator Source code: share/lib/hoc/subiter.hoc
		// Note that the dimensions of the value are not accounted for yet !!!
		// The only parameter in use (distance) is passed as dimensionless microns.
		
		// the inhomogeneous parameter comes from a source metric, such as tha following:
		enum{
			PATH_LENGTH_FROM_ROOT
		}metric;
		
		// flags for pre-processing the source metric through the context of this segment group, 
		// before passing it to the function that eventually determines an inhomogeneous distribution
		
		// from the metric that is evaluated over the group, should the minimum be subtracted from all samples so that the new minimum is zero?
		bool subtract_the_minimum;
		// after posibly subtracting the minimum, should all samples be divided by the maximum so that the new maximum is one (if the samples are non zero )
		bool divide_by_maximum;
		
		// The name of the variable within the variableParameter formulae that use it. Not to be confused with the id of the inhomogeneousParameter. 
		std::string variable;
	};
	
	struct SegmentGroup{
		// List of segments that the group eventually contains
		IdListRle list;
		
		// The group could be be a cable, in which case:
        //    This group contains an unbranched set of segments, and all of the segmentGroups marked with
        //    neuroLexId = sao864921383 form a non-overlapping set of all of the segments. 
        //    These segmentGroups correspond to the 'cables' of NeuroML v1.8.1.
		bool is_cable; // unbranched section of neurite, represents neuroLexId="sao864921383" attribute
		Int cable_nseg; // number of compartments to use for this cable, -1 if unset (default number: 1 compartment for the whole section). Represents numberInternalDivisions property
		
		// Note: In the current implementaion of the jNML exporter, it seems that the names of inhomogeneousParameters have global scope - they are not prefixed by morphology and segment group explicitly.
		// Maybe it is intended for the functional of an inhomogeneousParameter to be applied on other segment groups.
		// In that case, the inhomogeneousParaemters should either have a local scope inside the morphology or segmentGroup, 
		// or they should be accompanied by a mapping from a model-wide unique name for an inhomogeneousParameter to the morphology and segmentGroup which contains it.
		
		CollectionWithNames<InhomogeneousParameter> inhomogeneous_parameters;
		
		void addd(Int id){ list.Addd(id); }
		void add(const SegmentGroup &group){ list.Add(group.list); }
		
		void debug_print() const {
			printf("Internal structure:\n");
			list.debug_print();
			printf("\n");
		}
		SegmentGroup(){
			is_cable = false;
			cable_nseg = -1;			
		}
	};
	std::vector<SegmentGroup> segment_groups;
	NameIndexer segment_groups_by_name; // TODO combine into CollectionWithNames
	
	//Get internal position in the segments array, from NeuroML id
	Int lookupSegId(Int nml_id) const{
		return segments.getSequential(nml_id);
	}
	Int lookupNmlId(Int seq_id) const{
		return segments.getId(seq_id);
	}
	bool addSegment(const Segment & new_seg, Int new_id){
		return segments.add(new_seg, new_id);
	}
	
	// containing internal seg ID's to preserve sanity
	IdListRle getFullList() const{
		IdListRle ret;
		ret.Addd(0, segments.size());
		return ret;
	}
	
	bool add(const SegmentGroup &new_group, const char *new_name){
		if(segment_groups_by_name.count(new_name)) return false;
			
		segment_groups.push_back(new_group);
		segment_groups_by_name.insert(std::make_pair(new_name, segment_groups.size()-1));
		return true;
	}
	std::string Stringify_SegSeq_List(const IdListRle &seg_seq_list) const{
		//Now this is tricky because the internal id's are not the same as. the NeuroML id's.
		IdListRle nml_id_list;
		seg_seq_list.reduce(
			[&](Int internal_id){
				nml_id_list.Addd(lookupNmlId(internal_id));
			}
		);
		nml_id_list.Compact();
		return nml_id_list.Stringify();
	}
	std::string Stringify(const SegmentGroup &group) const{
		return Stringify_SegSeq_List(group.list);
	}
	
	void debug_print() const {
		//debug output
		for(size_t i = 0; i < segments.contents.size(); i++){
			
			printf( "Segment %ld: ", lookupNmlId(i) );
			segments.atSeq(i).debug_print();
			
			printf("\n");
		}
		for(std::pair<const char *, Int> keyval : segment_groups_by_name){
			printf("Segment group %s: %s\n", keyval.first, Stringify(segment_groups[keyval.second]).c_str() );
		}
	}
};

// Specifiers for biophysical parameters across parts of a physical cell
struct AcrossSegOrSegGroup{
	enum Type{
		NONE,
		SEGMENT, // just one segment
		GROUP // a group of segments
	} type;
	Int seqid; // internal sequence-type ID of segment, or group
	
	AcrossSegOrSegGroup(){type = NONE; seqid = -1;}
	void Segment(Int _i){type = SEGMENT; seqid = _i;}
	void Group(Int _i){type = GROUP; seqid = _i;}
	
	// Receives a functor that is callable with resolved segment seqid's as argument
	template<typename Functor>
	Functor reduce(const Morphology &morph, Functor doit) const {
		if(type == SEGMENT){
			doit(seqid);
		}
		else if(type == GROUP){
			morph.segment_groups[seqid].list.reduce(doit);
		}
		else{
			// should not happen, assume empty set
		}
		
		return doit;
	}
	
	IdListRle toList(const Morphology &morph) const {
		IdListRle ret;
		if(type == SEGMENT){
			ret.Addd(seqid);
		}
		else if(type == GROUP){
			return morph.segment_groups[seqid].list;
		}
		else{
			// should not happen, assume empty set
		}
		
		return ret;
	}
	
	void debug_print(const Morphology &morph) const ;
};

struct ValueAcrossSegOrSegGroup : public AcrossSegOrSegGroup{
	Real value;
	//Fill a segment count-sized array
	// NOTE: perhaps this should be removed, since segments do not exactly map to compartments. So why else would one apply the spec?
	template<typename SegmentPropertyArray>
	void apply( const Morphology &morph, SegmentPropertyArray &arr ) const {
		AcrossSegOrSegGroup::reduce( morph, [ &arr, this ] (Int seqid) {
			arr[seqid] = this->value;
		});
	}
};

struct SpeciesAcrossSegOrSegGroup : public AcrossSegOrSegGroup{
	Int species;
	Int concentrationModel;
	Real initialConcentration;
	Real initialExtConcentration;
	
	std::string Stringify(const CollectionWithNames<IonSpecies> &ion_species) const ;
};

// All ion channel density specifications, rolled into one
struct ChannelDistribution : public AcrossSegOrSegGroup{
	
	// TODO this should be an unique identifier!
	const char *name; 
	
	//These are standard for all channel distributions
	Int ion_species;
	
	Int ion_channel;
	
	struct InhomogeneousValue{
		Int parm; // of associated morphology and segment group. The distribution must refer to a specific segment group to pass parsing
		
		// this syntax tree does not need its symbols resolved, since there is only one possible symbol: the referenced inhomogeneous parameter. 
		// (Even that could be missing when the expression is a constant)
		TermTable value;
	};
	
	enum Type{
		POPULATION,
		FIXED,
		VSHIFT,
		NERNST,
		NERNST_CA2,
		GHK,
		GHK2,
	} type;
	
	//Only a few of these may be used, in their relevant types
	// union{
		struct Conductivity{
			enum Type{
				NONE,
				FIXED,
				NON_UNIFORM
			} type;
			// union{
				Real value; // if fixed
				InhomogeneousValue inho; // if non-uniform
			// };
		} conductivity; // for most
		Int number; // for integral population of channels
		Real permeability; // for GHK1
	// };
	Real erev; // for Fixed and Population
	Real vshift; //for vshift
	
	// TODO put in the set of ions that control this Nernst/GHK sort of Erev somewhere, when NeuroML gets to support this
};

struct ConcentrationModel{
	//Just two, nearly the same, for now
	enum Type{
		LEAKY,
		FIXED_FACTOR,
		COMPONENT
	} type;
	
	Int ion_species;
	Real restingConc; // steady state
	Real decayConstant; // tau of exponential decay
	Real shellThickness_or_rhoFactor;
	
	ComponentInstance component;
};

struct Q10Settings{
	
	enum Type{
		FIXED, // regardless of temperature
		FACTOR // = factor ^ ( (temperature - experimental) / 10 )
	} type;
	Real q10;
	Real experimentalTemp;
	
};

struct IonChannel{
		
	struct Rate{
		//Also Variable, Time, and whatever else it needs to be
		enum Type{
			EXPONENTIAL,
			EXPLINEAR,
			SIGMOID,
			FIXED, //the additional 'constant' factor
			COMPONENT
		} type;
		union{
			struct{
				Real rate;
				Real midpoint;
				Real scale;
				Real constant;
			} formula;
		};
		ComponentInstance component;
	};
	
	// NeuroML does not specify initial gate values; they are initialized to derived 'inf' value, through LEMS
	struct GateBase{
		Int instances; // positive integer gating power in the HH model
	};
	struct GateBaseDynamic : public GateBase{
		Q10Settings q10; // since there is no "none" option this value must be set always, default value is fixed q10 = 1 
	};
	
	struct GateHHRates : public GateBaseDynamic{
		Rate forwardRate;
		Rate reverseRate;
	};
	
	struct GateHHTauInf : public GateBaseDynamic{
		Rate timeCourse; 
		Rate steadyState;
	};
	/*
	struct GateHHRatesTau : public GateBaseDynamic{
		Rate forwardRate;
		Rate reverseRate;
		Rate tau; //and inf is computed from alpha, beta
	};
	
	struct GateHHRatesInf : public GateBaseDynamic{
		Rate forwardRate;
		Rate reverseRate;
		Rate inf; //and tau is computed from alpha, beta
	};*/
	
	struct GateHHRatesTauInf : public GateBaseDynamic{
		Rate forwardRate; // these are used nowhere xD
		Rate reverseRate; // except for indirect LEMS accesses, perhaps?
		Rate timeCourse; // tau, inf are
		Rate steadyState; // actually used
	};
	
	struct SubGateFractional : public GateHHTauInf{
		Real fraction_of_conductivity; // unitless
	};
	
	struct GateFractional : public GateBaseDynamic{
		std::vector<SubGateFractional> subgates;
	};
	
	struct GateInstantaneous : public GateBase{
		Rate steadyState;
	};
	
	struct GateKS : public GateBaseDynamic{
		
		NameIndexer state_names;
		std::vector<int> open_states;
		std::vector<int> closed_states;
		
		struct Transition_Base{
			Int from;
			Int to;
		};
		struct ForRevTransition : public Transition_Base{
			Rate forwardRate;
			Rate reverseRate;
		};
		
		struct TauInfTransition : public Transition_Base{
			Rate steadyState; // baseVoltageDepVariable Component
			Rate timeCourse; // baseVoltageDepTime Component
		};
		
		struct Transition{
			// refactor inheritance LATER
			enum Type{
				FORWARD_REVERSE,
				TAU_INF
			} type;
			//union{
				ForRevTransition forrev;
				TauInfTransition tauinf;
			//};
		
		};
		
		std::vector<Transition> transitions;
		
	};
	
	// Elide type to avoid special-casing copy-pasta in the code
	struct Gate{
		enum Type{
			NONE, //dummy value
			//PASSIVE, // ionChannelPassive
			
			RATES,
			TAUINF,
			RATESTAU,
			RATESINF,
			RATESTAUINF,
			
			INSTANTANEOUS,
			
			FRACTIONAL,
			KINETIC,
			
			COMPONENT
		}type;
		
		//union{
			GateHHRatesTauInf gaga;
			GateInstantaneous instantaneous;
			
			Int fractional;
			Int kinetic;
			
			ComponentInstance component;
		//};
	};
	
	struct ConductanceScaling{
		enum Type{
			NONE,
			Q10,
			COMPONENT
		} type;
		Q10Settings q10;
		ComponentInstance component;
		ConductanceScaling(){
			type = NONE;
		}
	};
	
	enum Type{
		NONE,
		NATIVE,
		COMPONENT
	};
	
	//all ion channels have these
	Type type;
	Int species; // Is -1 when non-specific
	Real conductance; // Is NaN when unspecified, is useful only for integer populations
	
	// and optionally this (allow only one to avoid madness until LATER)
	ConductanceScaling conductance_scaling;
	
	//and some of these, if a core NeuroML type
		CollectionWithNames<Gate> gates;
		// merge all the subgates of all composite ion channel gates to allow for union "polymorphism"
		std::vector<GateFractional> fractional_gates;
		std::vector<GateKS> kinetic_gates;
	
	// or it could be a Component
		ComponentInstance component;
};

// The biophysics specification of a physical cell
struct BiophysicalProperties{
	
	bool two_ca_pools;
	struct IntracellularProperties{
	
		std::vector<ValueAcrossSegOrSegGroup> resistivity_specs;
		std::vector<SpeciesAcrossSegOrSegGroup> ion_species_specs;
	
		IntracellularProperties(){  }
	}intracellularProperties;
	
	struct ExtracellularProperties{
	
		std::vector<SpeciesAcrossSegOrSegGroup> ion_species_specs; // what TODO with these anyway?
	
	}extracellularProperties;
	
	struct MembraneProperties{
		
		std::vector<ChannelDistribution> channel_specs;
		std::vector<ValueAcrossSegOrSegGroup> threshold_specs;
		std::vector<ValueAcrossSegOrSegGroup> capacitance_specs;
		std::vector<ValueAcrossSegOrSegGroup> initvolt_specs;
		
	}membraneProperties;
	
	Int morph_id;
	
	// NOTE this may be the same for different objects, since they may be different interpretations of the same <biophysicalProperties> tag, on different Morphologies.
	// Is an empty string if absent.
	const char * name; 
	
	// a small workaround to keep the 'calcium' name reference alive until persistent name tables are implemented LATER
	Int Ca_species_seq, Ca2_species_seq;
	
	void debug_print(const Morphology &morph, const CollectionWithNames<IonSpecies> &ion_species) const ;
};

// A model of a physical cells with shape, physiology and entirely physical properties
struct PhysicalCell : public Standalone{
	
	Int morphology;
	Int biophysicalProperties;

};

// an abstract(artificial) point neuron.
// Spike sources may also belong to this type, despite being inputs
struct ArtificialCell{
	enum Type{
		NONE = 0,
		
		IAF,
		IAF_REF,
		IAF_TAU,
		IAF_TAU_REF,
		IZH,
		IZH_2007,
		ADEX,
		FN,
		FN_1969,
		PINSKY_RINZEL_CA3,
		
		PYNN_IF_CURR_ALPHA,
		PYNN_IF_CURR_EXP,
		PYNN_IF_COND_ALPHA,
		PYNN_IF_COND_EXP,
		PYNN_EIF_COND_EXP_ISFA_ISTA,
		PYNN_EIF_COND_ALPHA_ISFA_ISTA,
		PYNN_HH_COND_EXP,
		
		COMPONENT, // non-core LEMS component
		SPIKE_SOURCE, // refer to input sources !
		
		MAX // sentinel value
	} type;
	
	
	Real C; // Capacitance factor, in whatever dimension
	
	
	Real thresh, reset; // for I&F cells
	
	Real leakConductance, leakReversal; // for leaky cells
	Real tau; // for time-constant I&F cells
	
	Real refract; // for refractory-period cells
	
	Real a,b,c,d, v0, vpeak, k; // for Izhikevich cells
	
	Real vt, delt; // for ADEX cells
	
	Real I, phi, w0; // for FN cells
	
	// for Pinsky-Rinzel CA3 cells (why is this a core component?)
	Real iSoma, iDend, gNmda, gAmpa, gc, gLs, gLd, gNa, gKdr, gCa, gKahp, gKC, eNa, eCa, eK, eL, qd0, pp, alphac, betac, cm;
	
	// for PyNN cells
	Real v_init, v_rest, v_spike, tau_m, tau_w, i_offset, tau_syn_E, tau_syn_I, e_rev_E, e_rev_I ;
	Real v_offset, e_rev_K, e_rev_Na, e_rev_leak, g_leak, gbar_K, gbar_Na;
	
	// or it could be a Component (LEMS equivalent may also be filled in here)
	ComponentInstance component;
	
	// or it could be a pure spike input source !
	Int spike_source_seq;
	
	bool GetCurrentInputAndDimension( const CollectionWithNames<ComponentType> &component_types, Dimension &dim_out ) const;
	bool GetVoltageExposureAndDimension( const CollectionWithNames<ComponentType> &component_types, Dimension &dim_out ) const;
	
	bool HasSpikeIn ( const CollectionWithNames<ComponentType> &component_types ) const;
	bool HasSpikeOut( const CollectionWithNames<ComponentType> &component_types, const CollectionWithNames<InputSource> &input_sources ) const;
	
	ArtificialCell(){ type = NONE; }
};

// could be physical, or artificial
struct CellType{
	enum Type{
		NONE = 0,
		PHYSICAL,
		ARTIFICIAL
	} type;
	
	//one of two applies
	PhysicalCell physical;
	ArtificialCell artificial;
	
	bool GetCurrentInputAndDimension( const CollectionWithNames<ComponentType> &component_types, Dimension &dim_out ) const;
	bool GetVoltageExposureAndDimension( const CollectionWithNames<ComponentType> &component_types, Dimension &dim_out ) const;
	
	bool HasSpikeIn ( const CollectionWithNames<ComponentType> &component_types ) const;
	bool HasSpikeOut( const CollectionWithNames<ComponentType> &component_types, const CollectionWithNames<InputSource> &input_sources ) const;
	
	CellType(){ type = NONE; }
};

// A NeuroML network, that's the object of simulation
struct Network{
	
	struct Population{
		struct Instance{
			// Int id;
			Real x,y,z;
		};
		
		Int component_cell; // cell type internal ID
		// Instance ID's are increasing but may be non-sequential !
		// This means NeuroML ID's must be mapped to position in the instance array !
		// All internal references are resolved serial positions in the instance list, to preserve sanity
		CollectionWithIds<Instance> instances;
		
	};
	
	struct Projection{
		
		struct Connection{
			// TODO abolish different syn types among instances if they were never meant to be !
			
			enum Type{
				SPIKING,    // made up of a directional spike connection, affecting postsynaptic cell only
				ELECTRICAL, // made up of a bidirectional voltage-based component, affecting both cells
				CONTINUOUS  // made up of two components, affecting pre and post cell in whatever way
			} type;
			
			Int   preCell; // instance #, in internal sequential form !
			Int   preSegment; // segment #, in internal sequential form !
			Real  preFractionAlong;
			Int  postCell; // instance #, in internal sequential form !
			Int  postSegment; // segment #, in internal sequential form !
			Real postFractionAlong;
			Real weight; // is NaN if missing
			Real delay;  // is NaN if missing
			union{
				Int synapse; // synaptic component for one-way chemical or two-way electrical synapses
				struct{
					Int preComponent; // synaptic components
					Int postComponent; // until LEMS components are supported
				} continuous;
			};
			Connection(){ weight = delay = NAN; }
		};
		
		Int presynapticPopulation;
		Int postsynapticPopulation;
		
		CollectionWithIds<Connection> connections; 
	};
	
	struct Input{
		// the component (out of input types)
		Int component_type;
		
		// and the target 
		Int population; // the target population
		Int cell_instance;
		Int segment;
		Real fractionAlong;
		// TODO reuse SegmentLocator?
		
		Real weight; // NaN if missing
		Input(){ weight = NAN; }
	};
	struct InputList{
		Int component;
		Int population;
		CollectionWithIds<Int> input_instances_seq;
	};
	
	// experimental extension
	// A TimeSeriesReader is a source of one or more time series, whose data is retrieved from various URL's,
	// presented to the model as a collection of data points uniformly structured as collections of named scalars, evolving over time.
	// the corresponding EventSetReader can be moved to the <Simulation> unlike this one, because waveforms need to be accessible in the LemsPath space, for VariableReferences and such.
	struct TimeSeriesReader{
		// each element in the time series is a multidimensional vector, with sub-elements of different dimensions
		struct InputColumn{
			Dimension dimension;
			LemsUnit units;
		};
		// more of a "url reference" in that relative paths are also accepted
		std::string source_url, data_format; // NB: to be parsed by the backend, for now, standardize LATER
		CollectionWithNames<InputColumn> columns; // set of properties, same for each element
		Int instances; // number of elements of the time series
	};
	
	
	Real temperature;
	
	CollectionWithNames<Population> populations;
	CollectionWithNames<Projection> projections;
	
	std::vector<Input> inputs; // all of them, whether in an inputList or sprinkled around as explicitinput, in order of appearance
	// TODO separate the explicitInput for a cleaner internal model
	// TODO add reverse index as well, for better debugging
	CollectionWithNames<InputList> input_lists;
	
	// experimental extension
	CollectionWithNames<TimeSeriesReader> data_readers;
};

// A NeuroML simulation, targeting a specific network
struct Simulation{
	
	// network is assumed to be the only one used in the simulation
	// TODO extract LemsCellLocator
	struct LemsSegmentLocator{
		// common for all, anything not applicable is set to -1
		// all internal id's are sequential, to preserve sanity
		Int population;
		Int cell_instance;
		
		Int segment_seq; // XXX explain what happens to this if it's a synapse or sth else
		// TODO move this to segment-based -- or not?
		Real fractionAlong; // Note that fractionAlong is officially missing ! Assume 0.5 (in the middle of the segment)
		// TODO add extension to the parser for "seg.fractionAlong" format
		
		// default invalid values, just in case
		LemsSegmentLocator(){
			population = -1;
			cell_instance = -1;
			segment_seq = -1;
			fractionAlong  = NAN;
		}
		LemsSegmentLocator(Int p, Int c, Int s, Real f):population(p), cell_instance(c), segment_seq(s), fractionAlong(f){}
		bool ok() const{
			return ( population >= 0 && cell_instance >= 0 && segment_seq >= 0 && 0 <= fractionAlong && fractionAlong <= 1 );
		}
	};
	
	// reference to inside the instance of a LEMS component
	struct LemsInstanceQuantityPath{
		Int namespace_thing_seq;
		// may support multiple levels of inheritance, etc LATER
	};
	
	struct MayBeLemsInstanceQuantity{
		LemsInstanceQuantityPath lems_quantity_path;
	};
	
	struct SynapticComponentQuantityPath : public MayBeLemsInstanceQuantity{
		struct Block : public MayBeLemsInstanceQuantity{
			enum Type{
				NONE,
				LEMS // only lemsified exist for now
			}type;	
		};
		struct Plasticity : public MayBeLemsInstanceQuantity{
			enum Type{
				NONE,
				LEMS // only lemsified exist for now
			}type;	
		};
		
		enum Type{
			NONE,
			NATIVE,
			BLOCK, // for such children
			PLASTICITY,
			LEMS // for a(n immediate) LEMS component
		}type;
		
		enum NativeEntry{
			GBASE,
			EREV,
			TAU, // later add more if needed, for now the rest are lemsified
			G,
			 // common for all, might be added later
			// WEIGHT,
			// DELAY
		}native_entry;
		
		Block block;
		Plasticity plasticity;
		
		// nothing else to add, for now
	};
	
	struct InputInstanceQuantityPath : public MayBeLemsInstanceQuantity{
		enum Type{
			NONE,
			NATIVE,
			SYNAPSE,
			LEMS // for a LEMS component
		}type;
		
		enum NativeEntry{
			AMPLITUDE,
			DURATION,
			DELAY
		}native_entry;
		
		SynapticComponentQuantityPath synapse;// for synapse containing inputs
		// nothing else to add, for now
	};
	
	struct LemsQuantityPath : public LemsSegmentLocator{
		
		enum Type{
			NONE,
			CELL, // artificial point neurons (or perhaps per-cell state variables?)
			SEGMENT, // only voltage is here, i guess - but also ca and ca2 concentration can be accessed from here!
			CHANNEL, // inside a channel(distribution), which may be a LEMS component as well
			ION_POOL, // may be a LEMS component as well
			SYNAPSE,  // may be a LEMS component as well
			INPUT, // an alternative way to access input instances
			NETWORK, // quite silly but networks can be (extended by) LEMS components, may be used LATER
			DATAREADER, // experimental extenbion
			MAX // sentinel
		};
		
		struct CellPath : public MayBeLemsInstanceQuantity{
			
			enum Type{
				NONE,
				INPUT, // for a cell that's actually an input
				// NATIVE, they all have lems representations for now
				LEMS // for a LEMS component
			}type;
			
			InputInstanceQuantityPath input; // in case
			// nothing else to add, for now
		};
		
		struct SegmentPath{
			//select the variable (some hardcoded options available)
			enum Type{
				NONE,
				VOLTAGE,
				CALCIUM_INTRA,
				CALCIUM2_INTRA,
				CALCIUM_EXTRA,
				CALCIUM2_EXTRA,
				// possible values: iChannels iSyn surfaceArea etc.
			}type;
			// nothing else to add
		};
		
		struct ChannelPath : public MayBeLemsInstanceQuantity{
			enum Type{
				NONE,
				I, // when channel is population based?
				G,
				G_DENSITY,	// when channel is density-based!
				I_DENSITY,
				Q, // a channel gate variable
				SUBQ, // a channel gate sub-variable
				LEMS
			}type;
			
			Int distribution_seq; // which cell-wide distribution spec? required
			Int gate_seq; // in case of particular gate
			Int subgate_seq; // in case of KS or composite gates
			
			// TODO add access to lems gates, gate rates and such!
		};
		
		struct IonPoolPath : public MayBeLemsInstanceQuantity{
			enum Type{
				NONE,
				CONCENTRATION_INTRA,
				CONCENTRATION_EXTRA,
				LEMS
			}type;
			
			Int distribution_seq; // which cell-wide distribution spec? required
			Int ion_type_seq; // TODO 
		};
		
		struct SynapsePath : public SynapticComponentQuantityPath{
			enum Location{
				PRE, // for general graded synapses
				POST // for chemical synapses and general graded synapses
				// NOTE: "both" would be ambiguous here... let's keep the paths specific
			}location;
			
			Int proj_seq; // projection in the network
			Int conn_seq; // connection in the list
			
			// just for convenience tracing down the path...? 
			// TODO canonize in Model, if that useful...?
			Int GetSynSeq(const Network::Projection::Connection &conn) const {
				if(conn.type == Network::Projection::Connection::CONTINUOUS){
					if ( location == Simulation::LemsQuantityPath::SynapsePath::Location::PRE ) return conn.continuous.preComponent;
					else return conn.continuous.postComponent;
				}
				else return conn.synapse;
			}
		};
		
		// LemsSegmentLocator is not used for these!
		struct InputPath : public InputInstanceQuantityPath{
			
			Int list_seq; // list in the network
			Int inst_seq; // instance in the list
			// LATER handle explicit connections if needed
		};
		
		struct DataReaderPath : public InputInstanceQuantityPath{
			Int read_seq; // reader in the network
			Int inst_seq; // instance in the group
			Int colu_seq; // property in the element
		};
		
		Type type;
		
		union{
			CellPath cell;
			SegmentPath segment;
			ChannelPath channel;
			IonPoolPath pool;
			SynapsePath synapse;
			InputPath input;
			DataReaderPath reader;
		};
		
		// in the wider sense that it refers to something attached to a specific cell
		// TODO refactor in the face of synapses which connect two cells ... ?
		// TODO remove, this concept makes sense only for the current version of Eden I guess ...? move to aux tooling of Model?
		// bool RefersToCell() const {
		// 	if(
		// 		type == CELL
		// 		|| type == SEGMENT || type == CHANNEL || type == ION_POOL
		// 		|| type == SYNAPSE || type == INPUT
		// 	) return true;
		// 	else return false;
		
		// }
		LemsQuantityPath(){
			type = NONE;
		}
	};
	
	
	struct LemsInstanceEventPath{
		enum Type{
			NONE,
			IN,
			OUT
		} type;
		Int event_port_seq;
		// may support multiple levels of inheritance, etc LATER
	};
	
	struct MayBeLemsInstanceEventPath{
		LemsInstanceEventPath lems_event_path;
	};
	
	struct InputInstanceEventPath : public MayBeLemsInstanceEventPath{
		enum Type{
			NONE,
			LEMS // for a LEMS component
		}type;
		
		// nothing else to add, for now
	};
	
	struct LemsEventPath : public LemsSegmentLocator{
		enum Type{
			NONE,
			CELL, // artificial point neurons
			SEGMENT,
			MAX // sentinel
		};
		struct Cell : public MayBeLemsInstanceEventPath{
			// applicable only to LEMS, only LEMS applies here -- at least for now?
			enum Type{
				NONE,
				INPUT, // for a cell that's actually an input
				LEMS
			}type;
			
			InputInstanceEventPath input; // in case
			// nothing else to add
		};
		struct Segment{
			enum Type{
				NONE,
				SPIKE, // 'spike' is always defined for a point neuron or physical compartment
				LEMS,
				MAX // sentinel
			};
			
			Type type;
		};
		
		Type type;
		
		union{
			Cell cell;
			Segment segment;
		};
		
		LemsEventPath(){
			type = NONE;
		}
	};
	
	struct LoggerBase{
		std::string fileName;
	};
	//for state variable trajectories
	struct DataWriter : public LoggerBase{
		struct OutputColumn{
			
			LemsQuantityPath quantity;
		};
		CollectionWithNames<OutputColumn> output_columns;
	};
	//for spikes or other discrete events
	struct EventWriter : public LoggerBase{
		struct EventSelection{
			
			LemsEventPath selection; // the component's port, 
		};
		enum{
			TIME_ID,
			ID_TIME
		} format; // as seen in https://github.com/LEMS/jLEMS/blob/development/src/main/java/org/lemsml/jlems/core/type/simulation/EventWriter.java#L18
		CollectionWithNames<EventSelection> outputs;
	};
	
	// extensions!
	struct CustomSetup{
		struct Statement{
			enum Type{
				NONE = 0,
				POPULATION,
				PROJECTION,
				INPUT_LIST
			};
			Type type;
			
			Int group_seq;
			std::vector<int> items_seq; // items of the cell population, or input list, or synapse list. Special case: empty means "all" in the order they were listed (ie seq).
			bool multi_mode; // are values different per item?
			
			Int seggroup_seq; // for physical cell entries
			
			// TODO IdListRle?
			std::vector<int> segments_seq; // for physical cell entries, *only* if a segment group is not listed; otherwise it is ignored.
			// special case: if seggroup_seq *and* segments_seq is empty, then the whole cell is set.
			
			bool cable_mode; // is a cable being sampled across its length?
			
			bool on_pre, on_post; // for synapse specifier
			
			// use this object for its parts really - not the whole of it, be careful! TODO explain
			Simulation::LemsQuantityPath path;
			
			std::vector<std::vector<Real> > real_data; // in engine units; one or many rows. If cable mode is set, the first row is fractionAlongCable sample points
			std::vector<std::vector<Simulation::LemsQuantityPath> > ref_data; // for 'reference' values.
			
			Statement(){
				group_seq = -1;
				multi_mode = false;
				
				seggroup_seq = -1;
				cable_mode = false;
				
				on_pre = on_post = false;
			}
			
		};
		std::vector<Statement> statements;	
	};
	// An EventSetReader is a set of event sources, streamed from various sorts of URL's.
	struct EventSetReader{
		struct EventMapping{
			Int source_port;
			LemsEventPath destination;
		};
		std::string source_url; // NB: to be parsed by the backend, for now, standardize LATER
		Int number_of_ports;
		std::vector<LemsQuantityPath> event_destinations;
		// TODO
	};
	// ---> fields
	
	Real length; // duration, that's how it's called in LEMS
	Real step; // timestep, likewise
	
	Int target_network;
	
	Int seed;
	bool seed_defined;
	
	// DataDisplay, Record, EventRecord not applicable yet
	CollectionWithNames<DataWriter> data_writers;
	CollectionWithNames<EventWriter> event_writers;
	
	// extensions!
	CustomSetup custom_init;
	
	CollectionWithNames<EventSetReader> event_readers; // FIXME
	
	
	Simulation(){
		seed_defined = false;
	}
};

// everything that can be imported from NeuroML files
struct Model{
	
	DimensionSet dimensions;
	CollectionWithNames<ComponentType> component_types;
	CollectionWithNames<ComponentInstance> component_instances;
	
	
	CollectionWithNames<Morphology> morphologies;
	std::vector<BiophysicalProperties> biophysics;
	
	CollectionWithNames<IonSpecies> ion_species;
	CollectionWithNames<ConcentrationModel> conc_models;
	
	CollectionWithNames<IonChannel> ion_channels;
	
	CollectionWithNames<CellType> cell_types;
	
	CollectionWithNames<SynapticComponent> synaptic_components;
	
	CollectionWithNames<InputSource> input_sources;
	
	CollectionWithNames<Network> networks;
	
	CollectionWithNames<Simulation> simulations;
	
	Int target_simulation;
	
	bool GetLemsQuantityPathType_FromLems(const Simulation::LemsInstanceQuantityPath &path, const ComponentInstance &compinst, ComponentType::NamespaceThing::Type &type, Dimension &dimension) const ;
	bool GetLemsQuantityPathType_SynapticComponent(const Simulation::SynapticComponentQuantityPath &path, const SynapticComponent &syn, ComponentType::NamespaceThing::Type &type, Dimension &dimension) const;
	bool GetLemsQuantityPathType_InputInstance(const Simulation::InputInstanceQuantityPath &path, const InputSource &input, ComponentType::NamespaceThing::Type &type, Dimension &dimension) const;
	bool GetLemsQuantityPathType(const Network &net, const Simulation::LemsQuantityPath &path, ComponentType::NamespaceThing::Type &type, Dimension &dimension) const;
	
	bool LemsQuantityPathToString(const ComponentInstance &inst, const Simulation::LemsInstanceQuantityPath &path, std::string &ret) const;
	bool LemsQuantityPathToString(const SynapticComponent &syn, const Simulation::SynapticComponentQuantityPath &path, std::string &ret) const;
	bool LemsQuantityPathToString(const InputSource &input, const Simulation::InputInstanceQuantityPath &path, std::string &ret) const;
	bool LemsQuantityPathToString(const ArtificialCell &cell, const Simulation::LemsQuantityPath::CellPath &path, std::string &ret) const;
	bool LemsQuantityPathToString(const Network &net, const Simulation::LemsQuantityPath &path, std::string &ret) const;
	
	bool ParseLemsSegmentLocator(const ILogProxy &log, const std::vector<std::string> tokens, const Network &net, Simulation::LemsSegmentLocator &path, Int &tokens_consumed) const;
	bool ParseLemsQuantityPathInComponent(const ILogProxy &log, const ComponentInstance &instance, const std::vector<std::string> &tokens, Simulation::LemsInstanceQuantityPath &lems_instance_qty_path, Int &tokens_consumed ) const;
	bool ParseLemsQuantityPath_SynapticComponent(const ILogProxy &log, const SynapticComponent &syn, const std::vector<std::string> &tokens, Simulation::SynapticComponentQuantityPath &path, Int &tokens_consumed ) const;
	bool ParseLemsQuantityPath_InputInstance(const ILogProxy &log, const InputSource &input, const std::vector<std::string> &tokens, Simulation::InputInstanceQuantityPath &path, Int &tokens_consumed ) const;
	bool ParseLemsQuantityPath_ArtificialCell(const ILogProxy &log, const ArtificialCell &cell,  const std::vector<std::string> &tokens, Simulation::LemsQuantityPath::CellPath &path_cell, Int &tokens_consumed) const;
	bool ParseLemsQuantityPath_CellProperty(const ILogProxy &log, const CellType &cell_type, const std::vector<std::string> &tokens, Simulation::LemsQuantityPath &path, Int &tokens_consumed ) const;
	bool ParseLemsQuantityPath(const ILogProxy &log, const char *qty_str, const Network &net, Simulation::LemsQuantityPath &path) const;
	
	Model(){ target_simulation = -1; }
};

//------------------> Parsed representations of NeuroML entities end

// Helper to keep load context, to preserve access to the model's strings after the model is loaded 
struct NmlImportContext_Holder{
	NmlImportContext_Holder(); NmlImportContext_Holder &operator=( NmlImportContext_Holder &&); ~NmlImportContext_Holder(); //; 
	struct Impl; std::unique_ptr<Impl> impl;
};

// Get the NeuroML model, use NmlImportContext_Holder to keep its strings loaded also
bool ReadNeuroML(const char *filename, Model &model, bool entire_simulation, bool verbose = false, FILE *info_log = stdout, FILE *error_log = stderr);
bool ReadNeuroML(const char *filename, Model &model, bool entire_simulation, NmlImportContext_Holder &import_context, bool verbose = false, FILE *info_log = stdout, FILE *error_log = stderr);

#endif
