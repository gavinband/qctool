
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
	void zlib_compress( char const* buffer, char const* const end, std::vector< char >* dest ) ;

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

	// Uncompress the given data, which should consist of a 64-bit unsigned integer in little-endian
	// format specifying the size in bytes of the uncompressed data, followed by the zlib-compressed data,
	// into the destination.  The destination must be large enough to fit the uncompressed data.
	template< typename T >
	void zlib_uncompress( std::vector< char > const& source, std::vector< T >* dest ) {
		#if HAVE_ZLIB
			uLongf const source_size = source.size() ;
			uLongf dest_size = dest->size() * sizeof( T ) ;
			int result = uncompress(
				reinterpret_cast< Bytef* >( &dest->operator[]( 0 ) ),
				&dest_size,
				reinterpret_cast< Bytef const* >( &source[0] ),
				source_size
			) ;
			assert( result == Z_OK ) ;
			assert( dest_size % sizeof( T ) == 0 ) ;
			dest->resize( dest_size / sizeof( T )) ;
		#else
			assert( 0 ) ; // no zlib support.
		#endif
	}
}

#endif
