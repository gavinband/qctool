#ifndef __GTOOL__FILEUTIL_HPP__
#define __GTOOL__FILEUTIL_HPP__

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <iterator>
#include <iostream>
#include <GToolException.hpp>

typedef std::auto_ptr< std::istream > INPUT_FILE_PTR ;
typedef std::auto_ptr< std::ostream > OUTPUT_FILE_PTR ;

struct FileException: public GToolException
{
	FileException( std::string const& msg )
		: GToolException( msg )
	{}
} ;

//
enum FileCompressionType {e_None = 0x1, e_Gzip = 0x2, e_Bzip2 = 0x4, e_FileCompressionMask = 0xf } ;
enum FileModeType { e_TextMode = 0x10, e_BinaryMode = 0x20, e_FileModeMask = 0xf0 } ;

FileCompressionType determine_file_compression( std::string const& filename ) ;
FileModeType determine_file_mode( std::string const& filename ) ;

// Return a stream representing a given input file, optionally with decompression
INPUT_FILE_PTR
open_file_for_input( std::string const& filename, int mode_flags ) ;

// Return a stream representing a given input file, attempting to auto-detect any compression used.
INPUT_FILE_PTR
open_file_for_input( std::string const& filename ) ;

// Return a stream representing a given output file, optionally with compression.
OUTPUT_FILE_PTR
open_file_for_output( std::string const& filename, int mode_flags ) ;

// Return a stream representing a given output file, guessing compression option from the filename.
OUTPUT_FILE_PTR
open_file_for_output( std::string const& filename ) ;

bool exists( std::string const& filename ) ;
bool is_regular( std::string const& filename ) ;

#endif

