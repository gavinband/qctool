
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include <zlib.h>
#include "genfile/zlib.hpp"

namespace genfile {
	void zlib_compress( char const* buffer, char const* const end, std::vector< char >* dest, std::size_t const offset ) {
		assert( dest != 0 ) ;
		uLongf const source_size = ( end - buffer ) ;
		uLongf compressed_size = compressBound( source_size ) ;
		dest->resize( compressed_size + offset ) ;
		int result = compress2(
			reinterpret_cast< Bytef* >( const_cast< char* >( &( dest->operator[](0) ) + offset ) ),
			&compressed_size,
			reinterpret_cast< Bytef const* >( buffer ),
			source_size,
			Z_BEST_COMPRESSION
		) ;
		assert( result == Z_OK ) ;
		dest->resize( compressed_size + offset ) ;
	}
}
