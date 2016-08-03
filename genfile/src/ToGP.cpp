
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
		
		// Generate a mapping from genotypes (encoded in a packed bit representation)
		// to indices in a specific ordering.  The ordering is colex order of the allele
		// count representation of the genotypes.  This order has the property that
		// genotypes carrying alleles k+1, ... always appear later in the order
		// than genotypes only carrying the first k alleles, for each k.
		// 
		// The genotype representation given by this function packs all possible genotypes
		// into a 16-bit word, and uses the minimum number of bits needed to encode genotypes
		// given the ploidy.  This number of bits is returned as the first member of the
		// return value, while the second member is the mapping from genotypes
		// to the index in the ordering.  The third member is the 'back' mapping
		// from index to encoded genotype.
		//
		// Example: for ploidy 3 we use 2 bits per allele.  The count of the 
		// first allele is omitted from the bit representation, and so the order
		// in vcf GT, allele count, and bit representation is
		// 0/0/0	(3,0,0,0)	...00 00 00 (11)
		// 0/0/1	(2,1,0,0)	...00 00 01 (10)
		// 0/1/1	(1,2,0,0)	...00 00 10 (01)
		// 1/1/1	(0,3,0,0)	...00 00 11 (00)
		// 0/0/2	(2,0,1,0)	...00 01 00 (10)
		// 0/1/2	(1,1,1,0)	...00 01 01 (01)
		// 1/1/2	(0,2,1,0)	...00 01 10 (00)
		// 1/2/2	(0,1,2,0)	...00 10 01 (00)
		// 2/2/2	(0,0,2,0)	...00 11 00 (00)
		// 2/2/3	(2,0,0,1)	...01 00 00 (10)
		// (etc.)
		//
		// A maximum of 16bits are used, so for a given ploidy the maximum possible number
		// of alleles in this representation is ⌊16/(bits per allele)⌋, where the bits per allele
		// is computed to be large enough to store the ploidy, i.e.
		// bits_per_allele = ⌈log2(ploidy)⌉.
		//
		Enumeration enumerate_unphased_genotypes(
			std::size_t const ploidy
		) {
			assert( ploidy < 128 ) ;
			if( ploidy == 0 ) {
				return Enumeration(
					std::make_pair( 0, std::numeric_limits< std::size_t >::max() ),
					std::make_pair( std::vector< uint16_t >( 1, 0 ), std::vector< uint16_t >( 1, 0 ) )
				) ;
			}
			// compute no. bits required to store values up to ploidy
			uint32_t bitsPerAllele = 0 ;
			for( std::size_t a = ploidy; a != 0; a >>= 1, ++bitsPerAllele ) ;
			std::size_t const numberOfAlleles = 16 / bitsPerAllele ;
			uint16_t const bitMask = uint16_t( 0xFFFF ) >> ( 16 - bitsPerAllele ) ;
			std::vector< uint16_t > limits( numberOfAlleles, ploidy ) ;
			bool finished = false ;
			std::vector< uint16_t > forwardMapping( 65536, 65535 ) ;
			std::vector< uint16_t > backMapping( 65536, 65535 ) ;
			uint16_t currentEncoded = 0 ;
			for( std::size_t index = 0; !finished; ++index ) {
				forwardMapping[ currentEncoded ] = index ;
				backMapping[ index ] = currentEncoded ;
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
			return Enumeration(
				std::make_pair( bitsPerAllele, numberOfAlleles ),
				std::make_pair( forwardMapping, backMapping )
			) ;
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
