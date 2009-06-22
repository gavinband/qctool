#ifndef __GTOOL__FILEUTIL_HPP__
#define __GTOOL__FILEUTIL_HPP__

#include <iostream>
#include <algorithm>
#include <string>
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

// Read a set of things from a file (optionally gzipped).
template< typename Set >
Set read_set_from_file( std::string filename, int mode_flags ) {
	Set aSet ;
	INPUT_FILE_PTR aStream( open_file_for_input( filename, mode_flags )) ;
	if( aStream.get() && aStream->good() ) {
		std::copy( std::istream_iterator< typename Set::value_type >( *aStream ),
		           std::istream_iterator< typename Set::value_type >(),
		           std::inserter( aSet, aSet.end() )) ;
	}

	return aSet ;
}


template< typename Set >
struct FromFileSet: Set
{
	typedef typename Set::value_type value_type ;
	
	FromFileSet( std::string filename ) {
		INPUT_FILE_PTR aStream( open_file_for_input( filename )) ;
	
		value_type value ;

		while( (*aStream) >> value ) {
			this->insert( value ) ;
		}
		
		if( aStream->bad() )
			throw FileException( "FromFileSet: After reading entries, stream was bad." ) ;
	}
} ;


#endif

