
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENBIN_REFERENCE_IMPLEMENTATION_IMPL_HPP
#define GENBIN_REFERENCE_IMPLEMENTATION_IMPL_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <stdint.h>
#include "../../config.hpp"
#ifdef HAVE_ZLIB
#include <sstream>
#include <zlib.h>
#endif
#include "genfile/endianness_utils.hpp"
#include "genfile/snp_data_utils.hpp"
#include "genfile/get_set.hpp"
#include "genfile/Error.hpp"
#include "genfile/zlib.hpp"

namespace genfile {
	namespace bgen {
		namespace impl {
			double get_probability_conversion_factor( uint32_t flags ) ;

			template< typename FloatType >
			uint16_t convert_to_integer_representation( FloatType number, FloatType factor ) {
				number *= factor ;
				if( number < 0 ) {
					std::cerr << "!! genfile::bgen::convert_to_integer_representation: probability " << number << " is -ve, clamping.\n" ;
					number = 0 ;
				}
				else if( number >= 65535.5 ) {
					std::cerr << "!! genfile::bgen::convert_to_integer_representation: probability " << number << " is too large, clamping.\n" ;
					number = 65535.0 ;
				}
				return static_cast< uint16_t > ( std::floor( number + 0.5 ) ) ;
			}

			template< typename FloatType >
			FloatType convert_from_integer_representation( uint16_t number, FloatType factor ) {
				FloatType result = number ;
				result /= factor ;
				return result ;
			}
			
			template<
				typename GenotypeProbabilitySetter
			>
			void read_uncompressed_snp_probability_data(
				char const* buffer,
				char const* const end,
				double const probability_conversion_factor,
				uint32_t number_of_samples,
				GenotypeProbabilitySetter set_genotype_probabilities
			) {
				for ( uint32_t i = 0 ; i < number_of_samples ; ++i ) {
					uint16_t AA = 0, AB = 0, BB = 0 ;
					assert( end >= buffer + 6 ) ;
					buffer = genfile::read_little_endian_integer( buffer, end, &AA ) ;
					buffer = genfile::read_little_endian_integer( buffer, end, &AB ) ;
					buffer = genfile::read_little_endian_integer( buffer, end, &BB ) ;

					set_genotype_probabilities(
						i,
						convert_from_integer_representation( AA, probability_conversion_factor ),
						convert_from_integer_representation( AB, probability_conversion_factor ),
						convert_from_integer_representation( BB, probability_conversion_factor )
					) ;
				}
			}
			

			template< typename GenotypeProbabilityGetter >
			void write_uncompressed_snp_probability_data(
				char* buffer,
				char* const end,
				uint32_t const flags,
				uint32_t number_of_samples,
				GenotypeProbabilityGetter get_AA_probability,
				GenotypeProbabilityGetter get_AB_probability,
				GenotypeProbabilityGetter get_BB_probability
			) {
				double const factor = impl::get_probability_conversion_factor( flags ) ;
				for ( uint32_t i = 0 ; i < number_of_samples ; ++i ) {
					uint16_t
						AA = convert_to_integer_representation( get_AA_probability( i ), factor ),
						AB = convert_to_integer_representation( get_AB_probability( i ), factor ),
						BB = convert_to_integer_representation( get_BB_probability( i ), factor ) ;
					assert( ( buffer + 6 ) <= end ) ;
					buffer = genfile::write_little_endian_integer( buffer, end, AA ) ;
					buffer = genfile::write_little_endian_integer( buffer, end, AB ) ;
					buffer = genfile::write_little_endian_integer( buffer, end, BB ) ;
				}
			}

			template< typename GenotypeProbabilityGetter >
			void write_compressed_snp_probability_data(
				std::ostream& aStream,
				uint32_t const flags,
				uint32_t number_of_samples,
				GenotypeProbabilityGetter get_AA_probability,
				GenotypeProbabilityGetter get_AB_probability,
				GenotypeProbabilityGetter get_BB_probability
			) {
				#if !HAVE_ZLIB
					assert(0) ; // zlib is required for compression support.
				#else
					// Construct a buffer into which we will compress.
					uLongf uncompressed_data_size = (6 * number_of_samples) ;
					uLongf buffer_size = 12 + (1.1 * uncompressed_data_size) ;		// calculated according to zlib manual.
					std::vector< Bytef > compression_buffer( buffer_size ) ;
					std::vector< Bytef > uncompressed_buffer( uncompressed_data_size ) ;
					// get probability data, uncompressed.
					write_uncompressed_snp_probability_data(
						reinterpret_cast< char* >( &uncompressed_buffer[0] ),
						reinterpret_cast< char* const >( &uncompressed_buffer[0] + uncompressed_data_size ),
						flags,
						number_of_samples,
						get_AA_probability, get_AB_probability, get_BB_probability
					) ;
					// compress it
					int result = compress( &compression_buffer[0], &buffer_size, &uncompressed_buffer[0], uncompressed_data_size ) ;
					assert( result == Z_OK ) ;
					// and write it (buffer_size is now the compressed length of the data).
					uint32_t the_buffer_size = buffer_size ;
					genfile::write_little_endian_integer( aStream, the_buffer_size ) ;
					aStream.write( reinterpret_cast< char const* >( &compression_buffer[0] ), buffer_size ) ;
				#endif
			}

			template< typename GenotypeProbabilityGetter >
			void write_uncompressed_snp_probability_data(
				std::ostream& aStream,
				uint32_t const flags,
				uint32_t const number_of_samples,
				GenotypeProbabilityGetter get_AA_probability,
				GenotypeProbabilityGetter get_AB_probability,
				GenotypeProbabilityGetter get_BB_probability
			) {
				std::size_t uncompressed_data_size = (6 * number_of_samples) ;
				std::vector< char > uncompressed_buffer( uncompressed_data_size ) ;
				// get probability data, uncompressed.
				write_uncompressed_snp_probability_data(
					&uncompressed_buffer[0],
					&uncompressed_buffer[0] + uncompressed_data_size,
					flags,
					number_of_samples,
					get_AA_probability, get_AB_probability, get_BB_probability
				) ;
				// compress it

				aStream.write( &uncompressed_buffer[0], uncompressed_data_size ) ;
			}
		}

	}
}

#endif
