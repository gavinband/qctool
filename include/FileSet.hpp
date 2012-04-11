
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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
Set read_set_from_file( std::string const& filename, int mode_flags = e_TextMode ) {
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
			if( value.find_first_not_of( "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_-+=!@£%$^&*()?/\\\";:~`[]{}" ) != std::string::npos ) {
				throw InvalidFileFormatError( filename ) ;				
			}
			this->insert( value ) ;
		}

		aStream->peek() ;  // flag eof if we're there.

		if( !aStream->eof() ) {
			throw InvalidFileFormatError( filename ) ;
		}

		if( aStream->bad() ) {
			throw InvalidFileFormatError( filename ) ;
		}
	}
} ;


#endif

