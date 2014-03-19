#include <cassert>
#include <Rcpp.h>
#include "zlib.h"
#define HAVE_ZLIB 1

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
				std::cerr << "Not ok!\n" ;
				throw std::exception() ;
			}
			if( dest_size % sizeof( T ) != 0 ) {
				std::cerr << "Not right size!\n" ;
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
	std::cerr << "result is " << std::dec << L << " x " << N << ".\n" ;
	
	std::vector< char > buffer ;
	std::vector< char > compressed_buffer ;
	for( int i = 0; i < L; ++i ) {
		RawVector const& data = rawData[i] ;
		compressed_buffer.assign( data.begin(), data.end() ) ;
		buffer.resize( N + 100 ) ;
		zlib_uncompress( &compressed_buffer[0], compressed_buffer.size(), &buffer ) ;
		if( i < 10 ) {
			std::cerr  << std::dec << "done, compressed size was " << data.size() << ", buffer size is " << buffer.size() << ", expected " << ( 10 + N ) << ".\n" ;
			std::cerr << "i = "  << std::dec << i << ", first few bytes of data are:\ncompressed: " ;
			for( std::size_t j = 0; j < data.size()
			; ++j ) { std::cerr << std::hex << int( data[j] ) ; }
			std::cerr << "...\nuncompressed: " ;
			for( std::size_t j = 0; j < buffer.size(); ++j ) { std::cerr << std::hex << int( buffer[j] ) << " " ; }
			std::cerr << "...\n" ;
		}
		assert( buffer.size() == N + 10 ) ;
		if( buffer[0] == 's' && buffer[1] == 0 && buffer[2] == 7 && std::string( buffer.begin() + 3, buffer.begin() + 10 ) == "bitpack" ) {
			for( std::size_t j = 0; j < N; ++j ) {
				result[i,j] = int( buffer[j+10] ) ;
				if( i < 10 & j < 20 ) {
					std::cerr << int( buffer[j+10] ) << " " ;
				}
			}
			if( i < 10 ) { std::cerr << "...\n" ; } 
		} else {
			// unsupported format.
			if( i < 10 ) {			std::cerr << "Unsupported format!\n" ; }
			throw std::exception() ;
		}
	}
	return result ;
}

