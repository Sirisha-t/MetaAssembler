#ifndef __FILEIO_H__
#define __FILEIO_H__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * File IO handler.
 */
namespace fio
{
	inline std::string getFileExtension( std::string filename )
	{
		return filename.substr( filename.find_last_of(".") + 1 );
	}

	/** Open file. */
	inline void openFile( std::fstream &file,
						  const char* name,
						  std::_Ios_Openmode mode )
	{
		file.open(name, mode);
		if ( !file) {
			std::cerr << "Can't open " << name << "\n";
			exit (EXIT_FAILURE);
		}
	}

	template<typename Type, typename Size>
	void write( Type *A, Size s, const char *filename )
	{
		FILE *pFile;
		pFile = fopen(filename, "wb");
		if ( fwrite(A, sizeof(Type), (size_t)s, pFile) != (size_t)s ) {
			fprintf(stderr, "Cannot write to `%s': ", filename);
			exit(EXIT_FAILURE);
		}
		fclose(pFile);
	}

	template<typename Type>
	size_t read( Type *&A, const char *filename, int extra )
	{
		FILE *pFile = fopen( filename, "rb" );
		long fsize;
		fseek(pFile, 0, SEEK_END);
		fsize = ftell(pFile);
		if ( fsize == 0 ) {
			std::cout << "Empty file:" << filename << "\n";
			exit(EXIT_FAILURE);
		}
		rewind(pFile);

		size_t nsize = (size_t)fsize/sizeof(Type);
		A = new Type[nsize+extra];
		size_t r = fread( A, sizeof(Type), nsize, pFile );
		if(r != nsize)  {
		  std::cout << "Reading error: " << filename << std::endl;
		}
		fclose(pFile);
		return nsize;
	}
}

#endif
