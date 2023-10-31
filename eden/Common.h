#ifndef COMMON_H
#define COMMON_H


#ifdef _WIN32

// Minimum supported version: Windows XP, aka 5.1
#define WINVER       0x0501
#define _WIN32_WINNT 0x0501

#define WIN32_LEAN_AND_MEAN // because all sorts of useful keywords like PURE are defined with the full headers
#define NOMINMAX
#include <windows.h>
// and more namespace pollution
#undef IN
#undef OUT
#undef INOUT
#undef CONST
#endif

#ifdef __APPLE__
// use these lines for non-OSX LATER, if you're targeting an iPad or something you know how to handle the specifics 
// #include <TargetConditionals.h>
// #if defined(TARGET_OS_OSX) && (TARGET_OS_OSX)
// ...
// #endif
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdint.h>
#include <float.h>
#include <stddef.h>
#include <inttypes.h>

#define _USE_MATH_DEFINES // maybe more elegant LATER, but since it works ...
#include <cmath> // apparently including <algorithm> undefines ::isfinite() function, on ICC

// for model representation and more
// TODO keep only what is really necessary
#include <vector>
#include <functional>
#include <algorithm>
#include <utility>
#include <string.h>
#include <string>

//for advanced whining
#include <stdarg.h>

//for time & memory usage measurement
#include <unistd.h> 


// Non-standard helper routines
#define pow10(powah) pow(10,(powah))
#define stricmp strcasecmp

//------------------> OS-independent utilities
// May be missing if not suported by platform, though

double TimevalDeltaSec(const timeval &start, const timeval &end);


// Memory measurements on Linux
// TODO support on more platforms !
#ifdef 	__linux__
//Memory consumption getters return 0 for not available
int64_t getCurrentResidentSetBytes();
int64_t getPeakResidentSetBytes();
int64_t getCurrentHeapBytes();
#endif

struct RunMetaData{
	double config_time_sec;
	double init_time_sec;
	double run_time_sec;
	double save_time_sec;
	
	int64_t peak_resident_memory_bytes; //0 for unknown
	int64_t end_resident_memory_bytes; //0 for unknown
	RunMetaData(){
		config_time_sec = NAN;
		init_time_sec = NAN;
		run_time_sec = NAN;
		save_time_sec = NAN;
		
		peak_resident_memory_bytes = 0;
		end_resident_memory_bytes = 0;
	}
};


// A very fast and chaotic RNG
// straight from Wikipedia
class XorShiftMul{
private:
		uint64_t state[1];
public:
	// NB if xorshiftmul's state becomes 0 it outputs zero forever :D
	// force it to not happen, with the highest bit set
	// In this case, certain combinations of shift factors ensure a complete cycle through all non-zero values.
	// It has been proven through linear algebra, somehow.
	XorShiftMul(uint64_t _seed){state[0] = _seed | (1LL << 63);}
	uint64_t Get(){
		uint64_t x = state[0];
		x ^= x >> 12; // a
		x ^= x << 25; // b
		x ^= x >> 27; // c
		state[0] = x;
		return x * 0x2545F4914F6CDD1D;
	}
};


// do not use default accuracy when converting numerics to alpha !
// explicitly specify what the alpha is used for
template< typename T, typename = typename std::enable_if< std::is_integral<T>::value >::type >
std::string accurate_string( T val ){
	return std::to_string( val );
}
// https://stackoverflow.com/questions/5985068/how-can-i-hide-defined-but-not-used-warnings-in-gcc
static inline std::string accurate_string( float val ){
	char tmps[100];
	sprintf(tmps, "%.9g", val);
	return tmps;
}
static inline std::string accurate_string( double val ){
	char tmps[100];
	sprintf(tmps, "%.17g", val);
	return tmps;
}

// Made static in hope they will not be inlined, wherever it makes sense
// but TODO check performance to justify it first!
//Uses strtol, overwrites errno.
//Returns whether conversion was successful.
static inline bool StrToL( const char *str, long &L, bool whole_string = true){
	long ret;
	char *pEnd;
	errno = 0;
	ret = strtol(str, &pEnd, 10);
	if( errno ) return false; //the standard way of handling strtol etc.
	if(whole_string){
		if(*pEnd) return false;
	}
	L = ret;
	return true;
}
static inline bool StrToF( const char *str, float &F ){
	float ret;
	char *pEnd;
	errno = 0;
	ret = strtof(str, &pEnd);
	if(errno) return false; //the standard way of handling strtol etc.
	while( *pEnd ){
		if(!isspace(*pEnd)) return false;
		pEnd++;
	}
	F = ret;
	return true;
}

// Split URL into scheme and auth+path, if scheme is present (otherwise it's a "URL reference")
bool GetUrlScheme(const std::string &url, std::string &scheme, std::string &auth_path);

// Tokenize a string, as with String.split() in string-capable languages
std::vector<std::string> string_split(const std::string& str, const std::string& delim);

bool GetLineColumnFromFile(const char *filename, const ptrdiff_t file_byte_offset, long long &line, long long &column);

//A more structured way to complain 
//Accepts varargs just because string manipulation is clunky in C/C++, even clunkier than varargs
void ReportErrorInFile_Base(FILE *error_log, const char *filename, const ptrdiff_t file_byte_offset, const char *format, va_list args);
void ReportErrorInFile(FILE *error_log, const char *filename, const ptrdiff_t file_byte_offset, const char *format, ...);

//------------------> Utilities end

//------------------> Windows specific util routines
#ifdef _WIN32
// TODO replace printf with a wrapper, on UTF16 targets
std::string DescribeErrorCode_Windows(DWORD error_code);
#endif
//------------------> end Windows specific util routines

#endif
