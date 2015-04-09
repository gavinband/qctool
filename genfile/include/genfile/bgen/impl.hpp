
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
#include "genfile/bgen/types.hpp"
#include "genfile/string_utils/hex.hpp"

// #define DEBUG_BGEN_FORMAT 3

namespace genfile {
	namespace bgen {
		struct Context {
			Context():
				number_of_samples(0),
				number_of_variants(0),
				magic( "bgen" ),
				free_data( "" ),
				flags(0)
			{}
				
			Context( Context const& other ):
				number_of_samples( other.number_of_samples ),
				number_of_variants( other.number_of_variants ),
				magic( other.magic ),
				free_data( other.free_data ),
				flags( other.flags )
			{}

			Context& operator=( Context const& other ) {
				number_of_samples = other.number_of_samples ;
				number_of_variants = other.number_of_variants ;
				magic = other.magic ;
				free_data = other.free_data ;
				flags = other.flags ;
				return *this ;
			}
			
			uint32_t header_size() const { return free_data.size() + 20 ; }
		public:	
			uint32_t number_of_samples ;
			uint32_t number_of_variants ;
			std::string magic ;
			std::string free_data ;
			uint32_t flags ;
		} ;
		
		namespace impl {
			uint32_t get_flags( std::string const& version ) ;
			
			double get_probability_conversion_factor( uint32_t flags ) ;

			template< typename FloatType >
			uint16_t convert_to_integer_representation( FloatType number, FloatType factor ) {
				number *= factor ;
				number = std::min( std::max( number, 0.0 ), 65535.0 ) ;
				return static_cast< uint16_t > ( std::floor( number + 0.5 ) ) ;
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
					Context const& context,
					GenotypeProbabilityGetter get_AA_probability,
					GenotypeProbabilityGetter get_AB_probability,
					GenotypeProbabilityGetter get_BB_probability
				) {
					double const factor = impl::get_probability_conversion_factor( context.flags ) ;
					for ( uint32_t i = 0 ; i < context.number_of_samples ; ++i ) {
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
				void round_probs_to_scaled_simplex( double* p, std::size_t* index, std::size_t const n, int const number_of_bits ) ;

				// Write data encoding n probabilities, given in probs, that sum to 1,
				// starting at the given offset in data.
				// Only t
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
					Context const& context,
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
					buffer = genfile::write_little_endian_integer( buffer, end, context.number_of_samples ) ;
					// Write ploidy
					uint16_t const numberOfAlleles = 2 ;
					uint8_t const ploidy = 2 ;
					buffer = genfile::write_little_endian_integer( buffer, end, numberOfAlleles ) ;
					buffer = genfile::write_little_endian_integer( buffer, end, ploidy ) ;
					buffer = genfile::write_little_endian_integer( buffer, end, ploidy ) ;
					char* ploidy_p = buffer ;
					buffer += context.number_of_samples ;
					buffer = genfile::write_little_endian_integer( buffer, end, uint8_t( 0 ) ) ;
					buffer = genfile::write_little_endian_integer( buffer, end, uint8_t( number_of_bits ) ) ;
					// We use an array of three doubles to compute the rounded probabilities.
					// We use a single 64-bit integer to marshall the data to be written.
					double v[3] ;
					std::size_t index[3] ;
					uint64_t data = 0 ;
					std::size_t offset = 0 ;
					for( std::size_t i = 0; i < context.number_of_samples; ++i ) {
						v[0] = get_AA_probability(i) ;
						v[1] = get_AB_probability(i) ;
						v[2] = get_BB_probability(i) ;
						bool missing = ( v[0] + v[1] + v[2] == 0 ) ;
						uint8_t ploidy = 2 | ( missing ? 0x80 : 0 ) ;
						*(ploidy_p++) = ploidy ;

						if( !missing ) {
							round_probs_to_scaled_simplex( &v[0], &index[0], 3, number_of_bits ) ;
							buffer = write_scaled_probs( &data, &offset, &v[0], 3, number_of_bits, buffer, end ) ;
	#if DEBUG_BGEN_FORMAT
							std::cerr << "genfile::bgen::impl::v12::write_uncompressed_snp_probability_data(): scaled probs are:"
								<< rounded_v[0] << ", " << rounded_v[1] << ", " << rounded_v[2] << ".\n" ;
							std::cerr << "genfile::bgen::impl::v12::write_uncompressed_snp_probability_data(): after write, data = "
								<< string_utils::to_hex(
									reinterpret_cast< unsigned char const* >( &data ),
									reinterpret_cast< unsigned char const* >( &data ) + 8
								) << ".\n" ;
	#endif
						}
					}
					// Get any leftover bytes.
					if( offset > 0 ) {
						int const nBytes = (offset+7)/8 ;
#if DEBUG_BGEN_FORMAT
						std::cerr << "genfile::bgen::impl::v12::write_uncompressed_snp_probability_data(): final offset = "
							<< offset << ", number of bits = " << number_of_bits << ", writing " << nBytes << " final bytes (space = " << (end-buffer) << ".\n" ;
#endif
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
				Context const& context,
				GenotypeProbabilityGetter get_AA_probability,
				GenotypeProbabilityGetter get_AB_probability,
				GenotypeProbabilityGetter get_BB_probability,
				int const number_of_bits = 16
			) {
				uint32_t const layout = context.flags & e_Layout ;
				if( layout == e_v11Layout ) {
					buffer = v11::write_uncompressed_snp_probability_data(
						buffer,
						bufferEnd,
						context,
						get_AA_probability, get_AB_probability, get_BB_probability
					) ;
					assert( buffer == bufferEnd ) ;
				} else if( layout == e_v12Layout ) {
					buffer = v12::write_uncompressed_snp_probability_data(
						buffer,
						bufferEnd,
						context,
						get_AA_probability, get_AB_probability, get_BB_probability,
						number_of_bits
					) ;
				} else {
					assert(0) ;
				}
				return buffer ;
			}
			
		}
	}
}

#endif
