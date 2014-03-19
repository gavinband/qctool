#include <cassert>
#include <Rcpp.h>
#include "zlib.h"
#define HAVE_ZLIB 1

// #define DEBUG_rcpp_uncompress_bitpack_haplotypes
using namespace Rcpp;

namespace {
	// Uncompress the given data, symmetric with zlib_compress.
	// The destination must be large enough to fit the uncompressed data,
	// and it will be resized to exactly fit the uncompressed data.
	template< typename T >
	void zlib_uncompress( char const* begin, uLongf const source_size, std::vector< T >* dest ) {
		#if HAVE_ZLIB
			uLongf dest_size = dest->size() * sizeof( T ) ;
			int result = uncompress(
				reinterpret_cast< Bytef* >( &dest->operator[]( 0 ) ),
				&dest_size,
				reinterpret_cast< Bytef const* >( begin ),
				source_size
			) ;
			if( result != Z_OK ) {
				throw std::exception() ;
			}
			if( dest_size % sizeof( T ) != 0 ) {
				throw std::exception() ;
			}
			dest->resize( dest_size / sizeof( T )) ;
		#else
			assert( 0 ) ; // no zlib support.
		#endif
	}
}

// [[Rcpp::export]]
IntegerMatrix rcpp_uncompress_bitpack_haplotypes( List rawData, int N ) {
	std::cerr << "sizeof( Rbyte ) = " << sizeof( Rbyte ) << ".\n" ;
	int const L = rawData.size() ;
	IntegerMatrix result( L, N ) ;

#if DEBUG_rcpp_uncompress_bitpack_haplotypes
	std::cerr << "result is " << std::dec << L << " x " << N << ".\n" ;
#endif
	
	std::vector< char > buffer ;
	std::vector< char > compressed_buffer ;
	for( int i = 0; i < L; ++i ) {
		RawVector const& data = rawData[i] ;
		compressed_buffer.assign( data.begin(), data.end() ) ;
		buffer.resize( N + 100 ) ;
		zlib_uncompress( &compressed_buffer[0], compressed_buffer.size(), &buffer ) ;
		if( buffer.size() != N + 10 ) {
			throw std::exception() ;
		}
		if( buffer[0] == 's' && buffer[1] == 0 && buffer[2] == 7 && std::string( buffer.begin() + 3, buffer.begin() + 10 ) == "bitpack" ) {
			for( std::size_t j = 0; j < N; ++j ) {
				result(i,j) = int( buffer[j+10] ) ;
			}
		} else {
			throw std::exception() ;
		}
	}
	return result ;
}

