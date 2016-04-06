
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <stdint.h>
#include <strstream>
#include "genfile/ToGP.hpp"

// #define DEBUG_TOGP 1
namespace genfile {
	namespace impl {
		// Generate all possible genotypes for a given ploidy and up to 4 alleles.
		// The order is the order such that genotypes carrying at least one allele k,...K
		// later in the order than those carrying no alleles k,....,K
		// This order is especially useful because (for the given ploidy) the genotypes at 
		// biallelic variants are ordered in the same way as those at a triallelic variant
		// but carrying none of the third allele - thus the orders are nested.
		//
		// Genotypes are stored as uint16_t with 4 bits per allele.  This gives
		// up to four alleles and up to ploidy of 15.
		std::vector< uint16_t > enumerate_unphased_genotypes(
			std::size_t ploidy
		) {
			std::size_t const max_number_of_alleles = 4 ;
#if DEBUG_TOGP
			std::cerr << "compute_genotype_ordering( " << ploidy << ", " << max_number_of_alleles << "):\n" ;
#endif
			// We will encode the genotypes in  uint16_t allowing 
			// 4 bits per allele, that is, up to four alleles and ploidy from 0 to 15.
			assert( max_number_of_alleles <= 4 ) ;
			assert( ploidy < 16 ) ;
			std::vector< uint16_t > limits( 4, ploidy ) ;
			// we assign 4 bits per allele, i.e. max ploidy is 15.
			uint16_t currentEncoded = 0 ;
			bool finished = false ;
			std::vector< uint16_t > result( 65536, 65535 ) ;
			for( std::size_t index = 0; !finished; ++index ) {
				result[ currentEncoded ] = index ;
#if DEBUG_TOGP
				std::cerr << "Stored genotype:" ;
				for( std::size_t i = 0; i < 4; ++i ) {
					std::cerr << " " << ((currentEncoded >> (4*i)) & 0xF) ;
				}
				std::cerr << ", limits:" ;
				for( std::size_t i = 0; i < 4; ++i ) {
					std::cerr << " " << limits[i] ;
				}
				std::cerr << "\n" ;
#endif
				std::size_t j = 0 ;
				for( ; j < 4; ++j ) {
					uint16_t value = (currentEncoded >> (j*4)) & 0xF ;
					if( value < limits[ j ] ) {
						currentEncoded += uint16_t( 1 ) << (j*4) ;
						for( std::size_t k = 0; k < j; ++k ) {
							assert( limits[k] > 0 ) ;
							--limits[k] ;
						}
						break ;
					} else {
						// Reset it to zero.
						// Note that to get here all lower-order counts must be zero.
						currentEncoded &= ( 0xFFFF << ((j+1)*4)) | (0xFFFF >> ((4-j)*4)) ;
						for( std::size_t k = 0; k < j; ++k ) {
							limits[k] += value ;
						}
					}
				}
				if( j == 4 ) {
					finished = true ;
				}
			}
			return result ;
		}

		std::string format_call( uint16_t call ) {
			std::ostringstream out ;
			for( std::size_t k = 0; k < 4; ++k ) {
				out << (k>0 ? "/": "") << ((call >> (4*k)) & 0xF) ;
			}
			return out.str() ;
		}

		std::vector< std::vector< uint16_t > > GTToGPUnphasedBase::m_tables ;
	}
}
