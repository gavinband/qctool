
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENBIN_REFERENCE_IMPLEMENTATION_IMPL_HPP
#define GENBIN_REFERENCE_IMPLEMENTATION_IMPL_HPP

#include <iostream>
#include <vector>
#include <map>
#include <cassert>
#include <cmath>
#include <cstdint>
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
#include "genfile/bgen/types.hpp"

#define DEBUG_BGEN_FORMAT 3

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
				return static_cast< uint16_t > ( std::round( number ) ) ;
			}

			namespace v11 {
				template< typename FloatType >
				FloatType convert_from_integer_representation( uint16_t number, FloatType factor ) {
					FloatType result = number ;
					result /= factor ;
					return result ;
				}
			
				template< typename GenotypeProbabilityGetter >
				char* write_uncompressed_snp_probability_data(
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
					return buffer ;
				}
			}
			
			namespace v12 {
				void round_probs_to_scaled_simplex( double* p, std::size_t const n, int const number_of_bits ) ;

				char* write_scaled_probs(
					uint64_t* data,
					std::size_t* offset,
					double const* probs,
					std::size_t const n,
					int const number_of_bits,
					char* buffer,
					char* const end
				) ;
					
				template< typename GenotypeProbabilityGetter >
				char* write_uncompressed_snp_probability_data(
					char* buffer,
					char* const end,
					uint32_t const flags,
					uint32_t number_of_samples,
					GenotypeProbabilityGetter get_AA_probability,
					GenotypeProbabilityGetter get_AB_probability,
					GenotypeProbabilityGetter get_BB_probability,
					int const number_of_bits
				) {
					assert( number_of_bits > 0 ) ;
					assert( number_of_bits <= 32 ) ;
#if DEBUG_BGEN_FORMAT
					std::cerr << "genfile::bgen::impl::v12::write_uncompressed_snp_probability_data(): number_of_bits = " << number_of_bits << ", buffer = " << reinterpret_cast< void* >( buffer ) << ", (end-buffer) = " << (end-buffer) << ".\n" ;
#endif
					buffer = genfile::write_little_endian_integer( buffer, end, number_of_samples ) ;
					// Write ploidy
					uint16_t const numberOfAlleles = 2 ;
					uint8_t const ploidy = 2 ;
					buffer = genfile::write_little_endian_integer( buffer, end, numberOfAlleles ) ;
					buffer = genfile::write_little_endian_integer( buffer, end, ploidy ) ;
					buffer = genfile::write_little_endian_integer( buffer, end, ploidy ) ;
					for( std::size_t i = 0; i < number_of_samples; ++i ) {
						buffer = genfile::write_little_endian_integer( buffer, end, ploidy ) ;
					}
					buffer = genfile::write_little_endian_integer( buffer, end, uint8_t( 0 ) ) ;
					buffer = genfile::write_little_endian_integer( buffer, end, uint8_t( number_of_bits ) ) ;
					// We use an array of three doubles to compute the rounded probabilities.
					// We use a single 64-bit integer to marshall the data to be written.
					double v[3] ;
					uint64_t data = 0 ;
					std::size_t offset = 0 ;
					for( std::size_t i = 0; i < number_of_samples; ++i ) {
						v[0] = get_AA_probability(i) ;
						v[1] = get_AB_probability(i) ;
						v[2] = get_BB_probability(i) ;
						round_probs_to_scaled_simplex( &v[0], 3, number_of_bits ) ;
						buffer = write_scaled_probs( &data, &offset, &v[0], 2, number_of_bits, buffer, end ) ;
						// Now write them
					}
					// Get any leftover bytes.
					if( offset > 0 ) {
						int const nBytes = (offset+7)/8 ;
						assert( (buffer+nBytes) <= end ) ;
						buffer = std::copy(
							reinterpret_cast< char const* >( &data ),
							reinterpret_cast< char const* >( &data ) + nBytes,
							buffer
						) ;
					}
					return buffer ;
				}
			}


			template< typename GenotypeProbabilityGetter >
			char* write_uncompressed_snp_probability_data(
				char* buffer,
				char* const bufferEnd,
				uint32_t const flags,
				uint32_t const number_of_samples,
				GenotypeProbabilityGetter get_AA_probability,
				GenotypeProbabilityGetter get_AB_probability,
				GenotypeProbabilityGetter get_BB_probability,
				int const number_of_bits = 16
			) {
				uint32_t const layout = flags & e_Layout ;
				if( layout == e_v11Layout ) {
					buffer = v11::write_uncompressed_snp_probability_data(
						buffer,
						bufferEnd,
						flags,
						number_of_samples,
						get_AA_probability, get_AB_probability, get_BB_probability
					) ;
					assert( buffer == bufferEnd ) ;
				} else if( layout == e_v12Layout ) {
					char* originalBuffer = buffer ;
					buffer += 4 ;
					buffer = v12::write_uncompressed_snp_probability_data(
						buffer,
						bufferEnd,
						flags,
						number_of_samples,
						get_AA_probability, get_AB_probability, get_BB_probability,
						number_of_bits
					) ;
					genfile::write_little_endian_integer( originalBuffer, originalBuffer+4, uint32_t( buffer - originalBuffer - 4 )) ;
				} else {
					assert(0) ;
				}
				return buffer ;
			}
			
			template< typename GenotypeProbabilityGetter >
			void write_uncompressed_snp_probability_data(
				std::ostream& aStream,
				uint32_t const flags,
				uint32_t const number_of_samples,
				GenotypeProbabilityGetter get_AA_probability,
				GenotypeProbabilityGetter get_AB_probability,
				GenotypeProbabilityGetter get_BB_probability,
				int const number_of_bits
			) {
				uint32_t const layout = flags & e_Layout ;
				// Construct a buffer into which we will compress.
				uLongf uncompressed_data_size = 
					( layout == e_v11Layout )
						? (6 * number_of_samples)
						: ( 14 + number_of_samples + 8 * ( ( number_of_samples * number_of_bits * 2 ) / 8 ) ) ;

				std::vector< char > uncompressed_buffer( uncompressed_data_size ) ;
				char* buffer = &uncompressed_buffer[0] ;
				char* end = write_uncompressed_snp_probability_data(
					buffer,
					buffer + uncompressed_data_size,
					flags,
					number_of_samples,
					get_AA_probability, get_AB_probability, get_BB_probability,
					number_of_bits
				) ;
				aStream.write( buffer, end - buffer ) ;
			}

			template< typename GenotypeProbabilityGetter >
			void write_compressed_snp_probability_data(
				std::ostream& aStream,
				uint32_t const flags,
				uint32_t number_of_samples,
				GenotypeProbabilityGetter get_AA_probability,
				GenotypeProbabilityGetter get_AB_probability,
				GenotypeProbabilityGetter get_BB_probability,
				int const number_of_bits = 16
			) {
				#if !HAVE_ZLIB
					assert(0) ; // zlib is required for compression support.
				#else
					uint32_t const layout = flags & e_Layout ;
					// Construct a buffer into which we will compress.
					uLongf uncompressed_data_size = 
						( layout == e_v11Layout )
							? (6 * number_of_samples)
							: ( 4 + 8 + number_of_samples + 8 * ( ( number_of_samples * number_of_bits * 2 ) / 8 ) ) ;
						
					uLongf buffer_size = 12 + (1.1 * uncompressed_data_size) ;		// calculated according to zlib manual.
					std::vector< Bytef > compression_buffer( buffer_size ) ;
					std::vector< Bytef > uncompressed_buffer( uncompressed_data_size ) ;
					// get probability data, uncompressed.
					write_uncompressed_snp_probability_data(
						reinterpret_cast< char* >( &uncompressed_buffer[0] ),
						reinterpret_cast< char* const >( &uncompressed_buffer[0] + uncompressed_data_size ),
						flags,
						number_of_samples,
						get_AA_probability, get_AB_probability, get_BB_probability,
						number_of_bits
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
		}
	}
}

#endif
