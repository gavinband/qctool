
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
		// Genotypes are stored as uint16_t with the given number of bits per allele.
		std::pair< uint32_t, std::vector< uint16_t > > enumerate_unphased_genotypes(
			std::size_t const ploidy
		) {
			assert( ploidy < 128 ) ;
			if( ploidy == 0 ) {
				return std::make_pair(
					0,
					std::vector< uint16_t >( 1, 0 )
				) ;
			}
			// compute no. bits required to store values up to ploidy
			uint32_t bitsPerAllele = 0 ;
			for( std::size_t a = ploidy; a != 0; a >>= 1, ++bitsPerAllele ) ;
			std::size_t const numberOfAlleles = 16 / bitsPerAllele ;
			uint16_t const bitMask = uint16_t( 0xFFFF ) >> ( 16 - bitsPerAllele ) ;
			std::vector< uint16_t > limits( numberOfAlleles, ploidy ) ;
			bool finished = false ;
			std::vector< uint16_t > result( 65536, 65535 ) ;
			
			uint16_t currentEncoded = 0 ;
			for( std::size_t index = 0; !finished; ++index ) {
				result[ currentEncoded ] = index ;
#if DEBUG_TOGP
				std::cerr << "ploidy: " << ploidy << ", " ;
				std::cerr << "Stored genotype:" ;
				for( std::size_t i = 0; i < numberOfAlleles; ++i ) {
					std::cerr << " " << ((currentEncoded >> (bitsPerAllele*i)) & bitMask) ;
				}
				std::cerr << ", limits:" ;
				for( std::size_t i = 0; i < numberOfAlleles; ++i ) {
					std::cerr << " " << limits[i] ;
				}
				std::cerr << "\n" ;
#endif
				std::size_t j = 0 ;
				for( ; j < numberOfAlleles; ++j ) {
					uint16_t value = (currentEncoded >> (j*bitsPerAllele)) & bitMask ;
					if( value < limits[ j ] ) {
						// value has not reached its limit; increase it by one.
						// Accordingly, lower the limits for all lower values by one.
						currentEncoded += uint16_t( 1 ) << (j*bitsPerAllele) ;
						for( std::size_t k = 0; k < j; ++k ) {
							assert( limits[k] > 0 ) ;
							--limits[k] ;
						}
						break ;
					} else {
						// Value has reached its limit.
						// Reset it to zero and reset the limits for lower values accordingly.
						// Note to get here all the values for lower alleles must be zero.
						currentEncoded &= (( 0xFFFF << ((j+1)*bitsPerAllele)) | ~( 0xFFFF << (j*bitsPerAllele )));
						for( std::size_t k = 0; k < j; ++k ) {
							limits[k] += value ;
						}
					}
				}
				if( j == numberOfAlleles ) {
					finished = true ;
				}
			}
			return std::make_pair( bitsPerAllele, result ) ;
		}

		std::string format_call( uint16_t call, uint32_t const bitsPerAllele ) {
			std::ostringstream out ;
			uint16_t const bitMask = uint16_t( 0xFFFF ) >> ( 16 - bitsPerAllele ) ;
			std::size_t const numberOfAlleles = 16 / bitsPerAllele ;
			for( std::size_t k = 0; k < numberOfAlleles; ++k ) {
				out << (k>0 ? "/": "") << ((call >> (numberOfAlleles*k)) & bitMask) ;
			}
			return out.str() ;
		}

		std::vector< std::pair< uint32_t, std::vector< uint16_t > > > GTToGPUnphasedBase::m_tables ;
	}
}
