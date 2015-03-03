
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#ifndef BGEN_REFERENCE_IMPLEMENTATION_HPP
#define BGEN_REFERENCE_IMPLEMENTATION_HPP

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
#include "genfile/bgen/impl.hpp"
#include "genfile/bgen/types.hpp"
#include "genfile/VariantDataReader.hpp"

/*
* This file contains a reference implementation of the BGEN file format
* specification described at:
* http://www.well.ox.ac.uk/~gav/bgen_format/bgen_format.html
*
* To use this file you will also need "endianness_utils.hpp".
*
*/

namespace genfile {
	namespace bgen {
		// Read the offset from the start of the stream.
		void read_offset( std::istream& iStream, uint32_t* offset ) ;
		// Write an offset value to the stream.
		void write_offset( std::ostream& oStream, uint32_t const offset ) ;

		// Read a header block from the supplied stream,
		// filling the fields of the supplied context object.
		std::size_t read_header_block(
			std::istream& aStream,
			Context* context
		) ;

		// Write a bgen header block to the supplied stream,
		// taking data from the fields of the supplied context object.
		void write_header_block(
			std::ostream& aStream,
			Context const& context
		) ;

		// Read a sample identifier block from the given stream.
		// The setter object passed in must be a unary function or
		// function object that takes a string.  It will be called
		// with each sample identifier that is read, in the order they
		// occur in the bgen file.
		template< typename SampleSetter >
		std::size_t read_sample_identifier_block(
			std::istream& aStream,
			Context const& context,
			SampleSetter setter
		) {
			uint32_t block_size = 0 ;
			uint32_t number_of_samples = 0 ;
			uint16_t identifier_size ;
			std::string identifier ;
			std::size_t bytes_read = 0 ;

			genfile::read_little_endian_integer( aStream, &block_size ) ;
			genfile::read_little_endian_integer( aStream, &number_of_samples ) ;
			bytes_read += 8 ;
			assert( number_of_samples == context.number_of_samples ) ;

			for( uint32_t i = 0; i < number_of_samples; ++i ) {
				genfile::read_length_followed_by_data( aStream, &identifier_size, &identifier ) ;
				if( aStream ) {
					bytes_read += sizeof( identifier_size ) + identifier_size ;
					setter( identifier ) ;
				} else {
					throw BGenError() ;
				}
			}
			assert( bytes_read == block_size ) ;
			return bytes_read ;
		}
		
		// Write the sample identifiers contained in the
		// given vector to the stream.
		std::size_t write_sample_identifier_block(
			std::ostream& aStream,
			Context const& context,
			std::vector< std::string > const& sample_ids
		) ;

		// Attempt to read identifying data fields for the next variant in the file.
		// This function will return true if SNP data was successfully read or false if the initial
		// reads met an EOF.  It will throw an BGenError if only a subset of fields can be read before EOF.
		template<
			typename NumberOfAllelesSetter,
			typename AlleleSetter
		>
		bool read_snp_identifying_data(
			std::istream& aStream,
			Context const& context,
			std::string* SNPID,
			std::string* RSID,
			std::string* chromosome,
			uint32_t* SNP_position,
			NumberOfAllelesSetter set_number_of_alleles,
			AlleleSetter set_allele
		) {
			uint16_t SNPID_size = 0;
			uint16_t RSID_size = 0;
            uint16_t numberOfAlleles = 0 ;
			uint16_t chromosome_size = 0 ;
			uint32_t allele_size = 0;
			std::string allele ;
			uint32_t const layout = context.flags & e_Layout ;
            
			if( layout == e_v11Layout || layout == e_v10Layout ) {
				uint32_t number_of_samples ;
				read_little_endian_integer( aStream, &number_of_samples ) ;
				if( !aStream ) {
					return false ;
				}
				if( number_of_samples != context.number_of_samples ) {
					throw BGenError() ;
				}
			}
			std::string chromosome_string ;
			read_length_followed_by_data( aStream, &SNPID_size, SNPID ) ;
			if( layout == e_v12Layout && !aStream ) {
				return false ;
			}
			read_length_followed_by_data( aStream, &RSID_size, RSID ) ;
			read_length_followed_by_data( aStream, &chromosome_size, chromosome ) ;
			read_little_endian_integer( aStream, SNP_position ) ;
			if( layout == e_v12Layout ) {
				read_little_endian_integer( aStream, &numberOfAlleles ) ;
			} else {
				numberOfAlleles = 2 ;
			}
			set_number_of_alleles( numberOfAlleles ) ;
			for( uint16_t i = 0; i < numberOfAlleles; ++i ) {
				read_length_followed_by_data( aStream, &allele_size, &allele ) ;
				set_allele( i, allele ) ;
			}
			if( !aStream ) {
				throw BGenError() ;
			}
			return true ;
		}

		// Read identifying data fields for the next variant in the file, assuming 2 alleles.
		// This function forwards to the generic multi-allele version, above, and will throw
		// a BGenError if the number of alleles is different than 2.
		bool read_snp_identifying_data(
			std::istream& aStream,
			Context const& context,
			std::string* SNPID,
			std::string* RSID,
			std::string* chromosome,
			uint32_t* SNP_position,
			std::string* first_allele,
			std::string* second_allele
		) ;
			
		// Write identifying data fields for the given variant.
		void write_snp_identifying_data(
			std::ostream& aStream,
			Context const& context,
			unsigned char max_id_size,
			std::string SNPID,
			std::string RSID,
			unsigned char chromosome,
			uint32_t SNP_position,
			std::string first_allele,
			std::string second_allele
		) ;

		// Read the raw (compressed or uncompressed) data for a variant from the stream
		void read_raw_probability_data(
			std::istream& aStream,
			Context const& context,
			std::vector< char >* buffer1
		) ;

		// uncompress or swap probability data from buffer1
		// into buffer2.  This may destroy the contents of buffer1.
		void uncompress_raw_probability_data(
			Context const& context,
			std::vector< char > const& buffer1,
			std::vector< char >* buffer2
		) ;

		// Skip over probability data for the current variant
		void ignore_snp_probability_data(
			std::istream& aStream,
			Context const& context
		) ;
		
		void parse_probability_data(
			char const* buffer,
			char const* const end,
			Context const& context,
			VariantDataReader::PerSampleSetter& setter
		) ;

		// Utility function which wraps the above steps into one.
		// It takes these steps:
		// 1: reads the length of the data
		// 2: reads that many bytes from the stream into buffer1
		// 3: if data is compressed, reads the uncompressed size of the data.
		// 4: uncompressed the remaining data into buffer2
		// 5: calls parse_probability_data to parse it, returning values to the setter.
		// The buffers will be resizes to fit data as needed.
		void read_and_parse_probability_data(
			std::istream& aStream,
			Context const& context,
			VariantDataReader::PerSampleSetter& setter,
			std::vector< char >* buffer1,
			std::vector< char >* buffer2
		) ;

		template< typename GenotypeProbabilityGetter >
		void write_snp_probability_data(
			std::ostream& aStream,
			Context const& context,
			GenotypeProbabilityGetter get_AA_probability,
			GenotypeProbabilityGetter get_AB_probability,
			GenotypeProbabilityGetter get_BB_probability,
			int const number_of_bits,
			std::vector< char >* buffer,
			std::vector< char >* compression_buffer
		) {
			uint32_t const layout = context.flags & e_Layout ;
			// Write the data to a buffer, which we then compress and write to the stream
			uLongf uncompressed_data_size =
				( layout == e_v11Layout )
					? (6 * context.number_of_samples)
					: ( 10 + context.number_of_samples + ((( context.number_of_samples * number_of_bits * 2 )+7) / 8) ) ;

			buffer->resize( uncompressed_data_size ) ;
			char* p = impl::write_uncompressed_snp_probability_data(
				&(*buffer)[0],
				&(*buffer)[0] + uncompressed_data_size,
				context,
				get_AA_probability, get_AB_probability, get_BB_probability,
				number_of_bits
			) ;
			assert( p = &(*buffer)[0] + uncompressed_data_size ) ;

			if( context.flags & e_CompressedSNPBlocks ) {
	#if HAVE_ZLIB
				uLongf compression_buffer_size = 12 + (1.1 * uncompressed_data_size) ;		// calculated according to zlib manual.
				compression_buffer->resize( compression_buffer_size ) ;
				int result = compress(
					reinterpret_cast< Bytef* >( &(*compression_buffer)[0] ), &compression_buffer_size,
					reinterpret_cast< Bytef* >( &(*buffer)[0] ), uncompressed_data_size
				) ;
				assert( result == Z_OK ) ;
				// write total payload size (compression_buffer_size is now the compressed length of the data).
				// Account for a 4-byte uncompressed data size if we are in layout 1.2
				if( layout == e_v12Layout ) {
					genfile::write_little_endian_integer( aStream, uint32_t( compression_buffer_size ) + 4 ) ;
					genfile::write_little_endian_integer( aStream, uint32_t( uncompressed_data_size )) ;
				} else {
					genfile::write_little_endian_integer( aStream, uint32_t( compression_buffer_size ) ) ;
				}

				// and write the data
				aStream.write( &(*compression_buffer)[0], compression_buffer_size ) ;
	#else
				assert(0) ;
	#endif
			} else {
				if( layout == e_v12Layout ) {
					genfile::write_little_endian_integer( aStream, uint32_t( uncompressed_data_size )) ;
				}
				aStream.write( &(*buffer)[0], uncompressed_data_size ) ;
			}
		}
	}
}

#endif
