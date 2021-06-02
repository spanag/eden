#include "Common.h"

std::vector<std::string> string_split(const std::string& str, const std::string& delim){
	// straight from StackOverflow https://stackoverflow.com/a/37454181
    std::vector<std::string> tokens;
    size_t prev = 0, pos = 0;
    do{
        pos = str.find(delim, prev);
        if (pos == std::string::npos) pos = str.length();
        std::string token = str.substr(prev, pos-prev);
        // if (!token.empty()) 
			tokens.push_back(token);
        prev = pos + delim.length();
    } while (pos < str.length() && prev < str.length());
    return tokens;
}

double TimevalDeltaSec(const timeval &start, const timeval &end){
	return (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 0.000001;
}

int64_t getCurrentResidentSetBytes(){
	int64_t ret = 0;
	FILE* fp = fopen( "/proc/self/statm", "r" );
	if(!fp) return 0;
	
	//second number should be resident size, in pages
	int64_t resident_pages;
	if ( fscanf( fp, "%*s %" SCNd64, &resident_pages ) == 1 ){
		//input ok
		ret = resident_pages * sysconf( _SC_PAGESIZE);
	}
	fclose( fp );
	return ret;
}


bool GetLineColumnFromFile(const char *filename, const ptrdiff_t file_byte_offset, long long &line, long long &column){
	FILE* f = fopen(filename, "rb");
	if (!f) return false;
	
	line = 1;
	ptrdiff_t offset = 0;
	ptrdiff_t last_line = 0;

	char buffer[1024];
	size_t size;

	while ((size = fread(buffer, 1, sizeof(buffer), f)) > 0)
	{
		for (long long i = 0; i < (long long)size; ++i){
			if (buffer[i] == '\n'){
				
				//result.push_back(offset + i);
				
				auto next_line = offset + i + 1;
				if(next_line > file_byte_offset) goto DONE;
				else{
					last_line = next_line;
					line++;
				}
			}
		}
		offset += size;
	}

	DONE:
	
	fclose(f);
	
	column = file_byte_offset - last_line + 1;
	return true;
}
int64_t getPeakResidentSetBytes(){
	struct rusage rusage = {0};
	getrusage( RUSAGE_SELF, &rusage );
	// available from Linux 2.6.32 onward! OSX has this value in bytes not kilobytes!
	return rusage.ru_maxrss * 1024LL;
}
int64_t getCurrentHeapBytes(){
	//system("cat /proc/self/smaps");
	int64_t ret = 0;
	//since Linux 2.6.14
	FILE* fp = fopen( "/proc/self/smaps", "r" );
	if(!fp) return 0;
	//search for the [heap] section, text processing in C yay!
	char *buf = NULL;
	size_t bufsize = 0;
	char in[85];
	while( !feof(fp) ){
		if(getline(&buf, &bufsize, fp) < 1) goto CLEANUP;
		
		if( sscanf(buf, "%*s %*s %*s %*s %*s %84s\n",in) && !strcmp(in,"[heap]") ){
			//printf("found heap section! %d\n",(int)feof(fp));
			
			while(!feof(fp)){
				//if(getline(&buf, &bufsize, fp) < 1) goto CLEANUP;
				
				if(fscanf(fp,"%84s",in) && !strcmp(in, "Referenced:") ){
					int64_t referenced_kilobytes;
					//printf("referenced!\n");
					if( fscanf(fp,"%" SCNd64, &referenced_kilobytes) ){
						ret = referenced_kilobytes*1024;
						goto CLEANUP;
					}
				}
				//else printf("%s ",in);
			}
		}//else printf("%s ",buf);
	}
	
	CLEANUP:
	fclose(fp);
	if(buf) free(buf);
	return ret;
}

void ReportErrorInFile_Base(FILE *error_log, const char *filename, const ptrdiff_t file_byte_offset, const char *format, va_list args){
	//file complaint header
	long long line = -1, column = -1;
	if(filename){
		if(file_byte_offset >= 0){
			//Apparently converting byte offset to line:column is complicated (and involves parsing the file once more!)	
			GetLineColumnFromFile(filename, file_byte_offset, line, column);
		}
		
		if(column > 0){
			fprintf(error_log, "file %s, line %lld, column %lld: ", filename, line, column);
		}
		else if(file_byte_offset > 0){
			fprintf(error_log, "file %s, offset %s: ", filename, std::to_string(file_byte_offset).c_str() );
		}
		else fprintf(error_log, "file %s: ", filename);
	}
	else{
		//unknown file
		if(file_byte_offset >= 0) fprintf(error_log, "<unknown file>, offset %s: ", std::to_string(file_byte_offset).c_str());
		else fprintf(error_log, "<unknown file>: ");
	}
	//specific error description
	vfprintf(error_log, format, args);
	//and a newline, why not
	fprintf(error_log, "\n");
	// and perhaps show the line
	auto ShowTheLine = []( const char *filename, long long file_byte_offset, long long column, FILE *error_log ){
		if( filename && column > 0 ){
			
			FILE *fp = fopen( filename, "r" );
			if( !fp ) return false;
			
			fseek( fp, file_byte_offset - (column - 1), SEEK_SET );
			char *line = NULL;
			size_t len = 0;
			if( getline(&line, &len, fp) < 0 ) return false;
			fprintf( error_log, "%s", line);
			for(int i = 0; i < column - 1 ;i++){
				if( line[i] == '\t' ) fprintf( error_log, ">\t" ); // for simple tab-correctness, more sophisticated LATER if ever needed
				else fprintf( error_log, "-" );
			}
			fprintf( error_log, "^\n" ); // end of arrow
			
			return true;
		}
		else return false;
	};
	ShowTheLine( filename, file_byte_offset, column, error_log  );
	
}
void ReportErrorInFile(FILE *error_log, const char *filename, const ptrdiff_t file_byte_offset, const char *format, ...){
	va_list args;
	va_start(args, format);
	ReportErrorInFile_Base(error_log, filename, file_byte_offset, format, args);
	va_end (args);
}
