#ifndef __GTOOL__FILESET_HPP__
#define __GTOOL__FILESET_HPP__

#include <iostream>
#include <algorithm>
#include <string>
#include <iterator>
#include <iostream>
#include <GToolException.hpp>
#include <FileUtil.hpp>

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
		if( !(*aStream) ) {
			throw FileException( "FromFileSet: Error opening file \"" + filename + "\".  Must be a readable file." ) ;
		}

		value_type value ;

		while( (*aStream) >> value ) {
			this->insert( value ) ;
		}
		
		if( aStream->bad() ) {
			throw FileException( "FromFileSet: Error reading entries -- the file must be a whitespace-separated list." ) ;
		}
	}
} ;


#endif

