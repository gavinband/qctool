
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_ZLIB_HPP
#define GENFILE_ZLIB_HPP

#include <vector>
#include <stdint.h>
#include <cassert>
#include "../config.hpp"
#ifdef HAVE_ZLIB
#include <sstream>
#include <zlib.h>
#endif
#include "genfile/Error.hpp"

namespace genfile {

	// Compress the given data into the given destination buffer.  The destination will be resized
	// to fit the compressed data.  (Since the capacity of dest may be larger than its size,
	// to save memory you may need to copy the contents of dest elsewhere after calling
	// this function).
	//
	// If offset is nonzero, compressed data will be written starting at position [offset].
	// The first [offset] bytes will be untouched.
	void zlib_compress( char const* buffer, char const* const end, std::vector< char >* dest, std::size_t const offset = 0 ) ;

	// Compress the given data into the given destination buffer.  The destination will be resized
	// to fit the compressed data.  (Since the capacity of dest may be larger than its size,
	// to save memory you may need to copy the contents of dest elsewhere after calling
	// this function).
	template< typename T >
	void zlib_compress( std::vector< T > const& source, std::vector< char >* dest ) {
		char const* begin = reinterpret_cast< char const* >( &source[0] ) ;
		char const* const end = reinterpret_cast< char const* >( &source[0] + source.size() ) ;
		return zlib_compress( begin, end, dest ) ;
	}

	template< typename T >
	void zlib_uncompress( char const* begin, char const* const end, std::vector< T >* dest ) {
	#if HAVE_ZLIB
		uLongf const source_size = ( end - begin ) ;
		uLongf dest_size = dest->size() * sizeof( T ) ;
		int result = uncompress(
			reinterpret_cast< Bytef* >( &dest->operator[]( 0 ) ),
			&dest_size,
			reinterpret_cast< Bytef const* >( begin ),
			source_size
		) ;
		assert( result == Z_OK ) ;
		assert( dest_size % sizeof( T ) == 0 ) ;
		dest->resize( dest_size / sizeof( T )) ;
	#else
		assert( 0 ) ; // no zlib support.
	#endif
	}

	// Uncompress the given data, symmetric with zlib_compress.
	// The destination must be large enough to fit the uncompressed data,
	// and it will be resized to exactly fit the uncompressed data.
	template< typename T >
	void zlib_uncompress( std::vector< char > const& source, std::vector< T >* dest ) {
		char const* begin = &source[0] ;
		char const* const end = &source[0] + source.size() ;
		zlib_uncompress( begin, end, dest ) ;
	}
}

#endif
