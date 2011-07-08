#ifndef GENFILE_ZLIB_HPP
#define GENFILE_ZLIB_HPP

#include <vector>
#include "../config.hpp"
#ifdef HAVE_ZLIB
#include <sstream>
#include <zlib.h>
#endif

namespace genfile {

	// Compress the given data into the given destination buffer, which will be resized
	// to fit the compressed data.  (Since the capacity of dest may be larger than its size,
	// to save memory it is advisable to copy the contents of dest elsewhere after calling
	// this function).
	template< typename T >
	void zlib_compress( std::vector< T > const& source, std::vector< char >* dest ) {
		assert( dest != 0 ) ;
		#if HAVE_ZLIB
			uLongf const source_size = source.size() * sizeof( T ) ;
			dest->resize( compressBound( source_size )) ;
			uLongf compressed_size = dest->size() ;
			int result = compress2(
				reinterpret_cast< Bytef* >( &dest->operator[]( 0 )),
				&compressed_size,
				reinterpret_cast< Bytef const* >( &source[0] ),
				source_size,
				Z_DEFAULT_COMPRESSION
			) ;
			assert( result == Z_OK ) ;
			dest->resize( compressed_size ) ;
		#else
			assert( 0 ) ; // no zlib support.
		#endif
	}

	// Uncompress the given data into the destination.
	// The destination buffer must be large enough to hold the uncompressed result.  It will be resized 
	// to fit the uncompressed data.
	template< typename T >
	void zlib_uncompress( std::vector< char > const& source, std::vector< T >* dest ) {
		#if HAVE_ZLIB
			uLongf const source_size = source.size() ;
			uLongf dest_size = dest->size() * sizeof( T ) ;
			int result = uncompress( reinterpret_cast< Bytef* >( &dest->operator[]( 0 ) ), &dest_size, reinterpret_cast< Bytef const* >( &source[0] ), source_size ) ;
			assert( result == Z_OK ) ;
			assert( dest_size % sizeof( T ) == 0 ) ;
			dest->resize( dest_size / sizeof( T )) ;
		#else
			assert( 0 ) ; // no zlib support.
		#endif
	}
}

#endif
