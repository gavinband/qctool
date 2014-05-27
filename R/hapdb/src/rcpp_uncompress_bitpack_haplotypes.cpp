#include <cassert>
#include <Rcpp.h>
#include "zlib.h"
#define HAVE_ZLIB 1

// #define DEBUG_rcpp_uncompress_bitpack_haplotypes 1
using namespace Rcpp;

namespace {
	// Uncompress the given data, symmetric with zlib_compress.
	// The destination must be large enough to fit the uncompressed data,
	// and it will be resized to exactly fit the uncompressed data.
	template< typename T >
	void zlib_uncompress( Rbyte const* begin, uLongf const source_size, std::vector< T >* dest ) {
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
	int const L = rawData.size() ;
	int const NN = N*2 ; // two haps per sample
	IntegerMatrix result( L, NN ) ;

#if DEBUG_rcpp_uncompress_bitpack_haplotypes
	std::cerr << "result is " << std::dec << result.nrow() << " x " << result.ncol() << ".\n" ;
#endif
	
	std::vector< char > buffer ;
	for( int i = 0; i < L; ++i ) {
		RawVector const& data = rawData[i] ;
		buffer.resize( NN + 10 ) ;
		zlib_uncompress( data.begin(), data.size(), &buffer ) ;
		if( buffer.size() != NN + 10 ) {
#if DEBUG_rcpp_uncompress_bitpack_haplotypes
			std::cerr << "Expected buffer size " << (NN+10) << " but got " << buffer.size() << ".\n" ;
#endif
			throw std::exception() ;
		}
		if( buffer[0] == 's' && buffer[1] == 0 && buffer[2] == 7 && std::string( buffer.begin() + 3, buffer.begin() + 10 ) == "bitpack" ) {
			for( std::size_t j = 0; j < NN; ++j ) {
				char const v = buffer[j+10] ;
				if( v == -1 ) {
					result(i,j) = NA_INTEGER ;
				} else {
					result(i,j) = int( v ) ;
				}
			}
		} else {
			throw std::exception() ;
		}
	}
	return result ;
}

// [[Rcpp::export]]
List rcpp_uncompress_floatarray_genotypes( List rawData, int N, IntegerVector const& chosen_samples, bool compute_probabilities, bool compute_dosage ) {
#if DEBUG_rcpp_uncompress_bitpack_haplotypes
	std::cerr << "rcpp_uncompress_floatarray_genotypes()\n" << std::flush ;
#endif
	int const L = rawData.size() ;
	
#if DEBUG_rcpp_uncompress_bitpack_haplotypes
	std::cerr << "result is " << std::dec << L << " x " << N << ".\n" ;
#endif

	NumericMatrix probabilities( 0, 0 ) ;
	NumericMatrix dosage( 0, 0 ) ;

	if( compute_probabilities ) {
		probabilities = NumericMatrix( L, chosen_samples.size()*3 ) ;
	}
	if( compute_dosage ){
		dosage = NumericMatrix( L, chosen_samples.size() ) ;
	}

#if DEBUG_rcpp_uncompress_bitpack_haplotypes
	std::cerr << "probabilities is " << std::dec << probabilities.nrow() << " x " << probabilities.ncol() << ".\n" ;
	std::cerr << "dosage is " << std::dec << dosage.nrow() << " x " << dosage.ncol() << ".\n" ;
#endif
	
	std::vector< unsigned char > buffer ;
	for( int i = 0; i < L; ++i ) {
		RawVector const& data = rawData[i] ;
		buffer.resize( (3*N)*4 + 21 ) ;
#if DEBUG_rcpp_uncompress_bitpack_haplotypes
	std::cerr << "data.size() == " << data.size() << ", uncompressed buffer size = " << buffer.size() << ".\n" ;
#endif
		zlib_uncompress( data.begin(), data.size(), &buffer ) ;
#if DEBUG_rcpp_uncompress_bitpack_haplotypes
	std::cerr << "data.size() == " << data.size() << ", uncompressed buffer size = " << buffer.size() << ".\n" ;
#endif
		if( buffer.size() != (3*N)*4 + 21 ) {
			throw std::exception() ;
		}
#if DEBUG_rcpp_uncompress_bitpack_haplotypes
	std::cerr << "data.size() == " << data.size() << ", uncompressed buffer size = " << buffer.size() << ".\n" ;
	std::cerr << "data: " << std::hex ;
	for( std::size_t index = 0; index < 21; ++index ) {
		std::cerr << int( buffer[index] ) << "(" << buffer[index] << ") " ;
	}
	std::cerr << std::dec << "\n" ;
#endif
		if( buffer[0] == 's' && buffer[1] == 0 && buffer[2] == 10 && std::string( buffer.begin() + 3, buffer.begin() + 13 ) == "floatarray" ) {
			uint32_t const N_samples = ( uint32_t( buffer[14] ) << 24 ) | ( uint32_t( buffer[15] ) << 16 ) | ( uint32_t( buffer[16] ) << 8 ) | ( uint32_t( buffer[17] ) << 0 ) ;
			uint16_t const N_entries_per_sample = ( uint16_t( buffer[19] ) << 8 ) | ( uint16_t( buffer[20] ) << 0 ) ;
#if DEBUG_rcpp_uncompress_bitpack_haplotypes
			std::cerr << "N_samples = " << N_samples << ", N_entries = " << N_entries_per_sample << ".\n" ;
#endif
			if( N_samples != N  ) {
				throw std::exception() ;
			}
			if( N_entries_per_sample != 3 ) {
				throw std::exception() ;
			}
			// ok, format looks good.  Let's read it:
			for( std::size_t ji = 0; ji < chosen_samples.size(); ++ji ) {
				std::size_t j = chosen_samples(ji) - 1 ; // Change R's 1-based indexing to 0-based
				if( compute_probabilities ) {
					for( std::size_t k = 0; k < 3; ++k ) {
						std::size_t const index = 21 + ((3*j)+k)*4 ;
						float const* value = reinterpret_cast< float const* >( (&buffer[0]) + index ) ;
						if( *value < 0 ) {
							probabilities(i,(3*j)+k) = NA_REAL ;
						} else {
							probabilities(i,(3*j)+k) = *value ;
						}
					}
				}
				
				if( compute_dosage ) {
					float const* g0 = reinterpret_cast< float const* >( (&buffer[0]) + 21 + ((3*j)+0)*4 ) ;
					float const* g1 = reinterpret_cast< float const* >( (&buffer[0]) + 21 + ((3*j)+1)*4 ) ;
					float const* g2 = reinterpret_cast< float const* >( (&buffer[0]) + 21 + ((3*j)+2)*4 ) ;
					if( *g0 < 0 | *g1 < 0 | *g2 < 0 ) {
						dosage(i,j) = NA_REAL ;
					} else {
						dosage(i,j) = *g1 + 2.0 * *g2 ;
					}
				}
			}
		} else {
			throw std::exception() ;
		}
	}
	return List::create( Rcpp::Named( "probabilities" ) = probabilities, Rcpp::Named( "dosage" ) = dosage ) ;
}

