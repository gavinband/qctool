#ifndef GENFILE_ZLIB_HPP
#define GENFILE_ZLIB_HPP

#include <vector>
#include <stdint.h>
#include "../config.hpp"
#ifdef HAVE_ZLIB
#include <sstream>
#include <zlib.h>
#endif
#include "genfile/Error.hpp"

namespace genfile {

	// Compress the given data into the given destination buffer, prepending with a 64-bit unsigned integer
	// in little-endian format recording the size in bytes of the uncompressed data.  The destination will be resized
	// to fit this integer plus the compressed data.  (Since the capacity of dest may be larger than its size,
	// to save memory it is advisable to copy the contents of dest elsewhere after calling
	// this function).
	template< typename T >
	void zlib_compress( std::vector< T > const& source, std::vector< char >* dest ) {
		assert( dest != 0 ) ;
		#if HAVE_ZLIB
			uLongf const source_size = source.size() * sizeof( T ) ;
			dest->resize( compressBound( source_size ) + 8 ) ;
			uLongf compressed_size = dest->size() ;
			int result = compress2(
				reinterpret_cast< Bytef* >( &dest->operator[]( 8 ) ),
				&compressed_size,
				reinterpret_cast< Bytef const* >( &source[0] ),
				source_size,
				Z_DEFAULT_COMPRESSION
			) ;
			assert( result == Z_OK ) ;
			assert( sizeof( uint64_t ) == 8 ) ;
			dest->resize( compressed_size + 8 ) ;
			// Write the compressed size as a 64-bit unsigned integer
			// in little-endian (low-order byte first) format.
			uint64_t compressed_size64 = compressed_size ;
			for( std::size_t byte = 0; byte < 8; ++byte ) {
				dest->operator[]( byte ) = ( compressed_size64 >> ( byte * 8 ) ) && 0xFF ;
			}
		#else
			assert( 0 ) ; // no zlib support.
		#endif
	}

	// Uncompress the given data, which should consist of a 64-bit unsigned integer in little-endian
	// format specifying the size in bytes of the uncompressed data, followed by the zlib-compressed data,
	// into the destination.  The destination will be resized to fit the uncompressed data.
	template< typename T >
	void zlib_uncompress( std::vector< char > const& source, std::vector< T >* dest ) {
		assert( source.size() >= 8 ) ;
		#if HAVE_ZLIB
			uint64_t uncompressed_size =
				( uint64_t( source[0] ) )
				+ ( uint64_t( source[1] ) << 8 )
				+ ( uint64_t( source[2] ) << 16 )
				+ ( uint64_t( source[3] ) << 24 )
				+ ( uint64_t( source[4] ) << 32 )
				+ ( uint64_t( source[5] ) << 40 )
				+ ( uint64_t( source[6] ) << 48 )
				+ ( uint64_t( source[7] ) << 56 ) ;
			;
			dest.resize( uncompressed_size / sizeof( T )) ;
			uLongf const source_size = source.size() - 8 ;
			uLongf dest_size = dest->size() * sizeof( T ) ;
			int result = uncompress(
				reinterpret_cast< Bytef* >( &dest->operator[]( 0 ) ),
				&dest_size,
				reinterpret_cast< Bytef const* >( &source[8] ),
				source_size
			) ;
			assert( result == Z_OK ) ;
			assert( dest_size % sizeof( T ) == 0 ) ;
			assert( uint64_t( dest_size ) == uncompressed_size ) ;
			dest->resize( dest_size / sizeof( T )) ;
		#else
			assert( 0 ) ; // no zlib support.
		#endif
	}
}

#endif
