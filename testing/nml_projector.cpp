
#include "Common.h"
#include "stdio.h"
#include "ctype.h"
// #include "thirdparty/pugixml-1.9/pugixml.hpp"

int main( int argc, char **argv ){
	
	if( argc != 6 ){
		fprintf( stderr, "usage: nml_projector <filename> <neurons> <density> <contact segment> <gap junction type>\n" );
		return 2;
	}
	
	const char *filename = NULL;
	int neurons = 0;
	float density = -1;
	int contact_segment = -1;
	const char *sGapJunction = NULL;
	
	filename = argv[1];
	// de-facto validation for this	
	
	const char *sNeurons = argv[2];
	if(!(
		sscanf( sNeurons, "%d", &neurons ) == 1
		&& neurons > 0
	)){
		fprintf( stderr, "error: neurons should be a positive integer, not %s\n", sNeurons );
	}
	
	const char *sDensity = argv[3];
	if(!(
		sscanf( sDensity, "%f", &density ) == 1
		&& 0 <= density && density <= 1
	)){
		fprintf( stderr, "error: density should be a number between 0 and 1, not %s\n", sDensity );
	}
	
	const char *sContactSegment = argv[4];
	if(!(
		sscanf( sContactSegment, "%d", &contact_segment ) == 1
		&& 0 <= contact_segment
	)){
		fprintf( stderr, "error: contact segment should be a non-negative integer, not %s\n", sContactSegment );
	}
	sGapJunction = argv[5];
	
	// no validation for this	
	
	auto get_bytes_from_file = []( const char *filename ){
		// release with free()
		FILE *file = NULL;
		char *filedata = NULL;
		size_t file_size = 0;
		size_t read_result = -1;
		
		bool ok = false;
		
		file = fopen( filename, "rb" );
		if(!file) goto CLEANUP;
		
		// get size
		if (fseek(file, 0, SEEK_END) != 0) goto CLEANUP;
		file_size = ftell(file);
		if (file_size < 0) goto CLEANUP;
		if (fseek(file, 0, SEEK_SET) != 0) goto CLEANUP;
		// load the file data
		filedata = (char *) malloc( file_size + 1 );
		if (!filedata) goto CLEANUP;
		// read the file into memory
		read_result = fread(filedata, sizeof(char), file_size, file);
		if (read_result != file_size){
			goto CLEANUP;
		}
		filedata[read_result] = '\0';
		
		ok = true;

		CLEANUP:
		if( !ok ){
			perror("loading file:");
		}
		if(file){
			fclose(file);
			file = NULL;
		}
			
		return filedata;
	};
	
	FILE *fout = NULL;
	char *filedata = NULL;
	
	filedata = get_bytes_from_file( filename );
	if( !filedata ){
		perror("Could not open file");
		return 1;
	}
	
	const char insert_before[] = "</network>";
	char *insert_position = strstr( filedata, insert_before );
	if( !insert_position ){
		printf( "%s not found in file", insert_before );
		return 3;
	}
	
	while( ( filedata < insert_position ) && isblank( *insert_position ) ) insert_position--;
	
	fout = fopen( filename, "wb" );
	if( !fout ){
		perror("writing file:");
		return 1;
	}
	
	fwrite( filedata, insert_position - filedata, 1, fout );
	
	const char sPre[] = "<electricalProjection id =\"gap_junctions\" presynapticPopulation=\"pop\" postsynapticPopulation=\"pop\">\n";
	const char sPost[] = "</electricalProjection>\n";
	
	fwrite( sPre, strlen(sPre), 1, fout );
	
	char tmps[10000];
	int conn_seq = 0;
	int random_seed = 0;
	
	const uint32_t sample_scale = (1 << 24); //how many bits wide should a random sample be?
	const uint32_t fraction_of_scale = uint32_t( density * sample_scale );
	for( int pre = 0; pre < neurons; pre++ ){
		XorShiftMul rng( ( (uint64_t)random_seed << 32) + pre);
		
		if( fraction_of_scale == 0 ) continue; // if completely empty, can skip populating
		for( int post = pre + 1; post < neurons; post++ ){
			
			uint32_t sample = rng.Get() % sample_scale;
			if(sample < fraction_of_scale){
				
				sprintf(tmps, "<electricalConnection id=\"%d\" preCell=\"%d\" postCell=\"%d\" preSegment=\"%d\" postSegment=\"%d\" synapse=\"%s\"/>\n", conn_seq, pre, post, contact_segment, contact_segment, sGapJunction);
				fwrite( tmps, strlen(tmps), 1, fout );
				
				conn_seq++;
			}
		}
	}
	
	fwrite( sPost, strlen(sPost), 1, fout );
	
	fwrite( insert_position, strlen(insert_position), 1, fout );
	fclose( fout );
	
	return 0;
}
