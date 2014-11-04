
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENBIN_REFERENCE_IMPLEMENTATION_HPP
#define GENBIN_REFERENCE_IMPLEMENTATION_HPP

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
		struct BGenError: public virtual std::exception {
			~BGenError() throw() {}
			char const* what() const throw() { return "BGenError" ; }
		} ;
		typedef ::uint32_t uint32_t ;
		typedef ::uint16_t uint16_t ;
		// v1.1 definitions were:
		// enum FlagsType { e_NoFlags = 0, e_CompressedSNPBlocks = 0x1, e_LongIds = 0x4 } ;
		// v1.2 definitions:
		enum FlagsType { e_NoFlags = 0, e_CompressedSNPBlocks = 0x1, e_Layout = 0x3C } ;
		enum Layout { e_v10Layout = 0x0, e_v11Layout = 0x4, e_v12Layout = 0x8 } ;

/*
		* Function: read_offset
		* Read the offset from the start of the stream.
*/
		template< typename OffsetType >
		void read_offset( std::istream& iStream, OffsetType* offset ) ;

/*
		* Function: read_offset
		* Read the offset value to the stream.
*/
		template< typename OffsetType >
		void write_offset( std::ostream& oStream, OffsetType offset ) ;

		/* Function: read_header_block
		* Read a header block from the supplied istream object, using
		* the supplied setter objects (or function pointers) to return
		* the data to the user.
*/
		template<
			typename HeaderSizeSetter,
			typename NumberOfSNPBlocksSetter,
			typename NumberOfSamplesSetter,
			typename SNPBlockSizeSetter,
			typename FlagsSetter,
			typename FreeDataSetter
		 >
		void read_header_block(
			std::istream& aStream,
			HeaderSizeSetter set_header_size,
			NumberOfSNPBlocksSetter set_number_of_snp_blocks,
			NumberOfSamplesSetter set_number_of_samples,
			FreeDataSetter set_free_data,
			FlagsSetter set_flags
		) ;

/*
		* Function: write_header_block()
		* Write a header block with the given information to the given ostream object.
*/
		void write_header_block(
			std::ostream& aStream,
			uint32_t number_of_snp_blocks,
			uint32_t number_of_samples,
			std::string const& free_data,
			uint32_t flags
		) ;

		/* Function: get_header_blocK-size()
		* Return the size in bytes of the header block with the given free data
*/
		std::size_t get_header_block_size(
			std::string const& free_data
		) ;

/*
		* Function: read_snp_block()
		* Read a snp block from the given istream object, returning the information
		* via the setter objects (or function pointers) passed as arguments.
		* These must be callable as:
		* - set_number_of_samples( integer )
		* - set_SNPID( string )
		* - set_RSID( string )
		* - set_SNP_position( integer )
		* - set_alleles( string, string )
		* - set_genotype_probabilities( double, double, double )
*/
		template<
			typename IntegerSetter,
			typename StringSetter,
			typename AlleleSetter,
			typename ChromosomeSetter,
			typename SNPPositionSetter,
			typename GenotypeProbabilitySetter
		 >
		void read_snp_block(
			std::istream& aStream,
			uint32_t const flags,
			IntegerSetter set_number_of_samples,
			StringSetter set_SNPID,
			StringSetter set_RSID,
			ChromosomeSetter set_chromosome,
			SNPPositionSetter set_SNP_position,
			AlleleSetter set_allele1,
			AlleleSetter set_allele2,
			GenotypeProbabilitySetter set_genotype_probabilities
		) ;

/*
		* Function: write_snp_block()
		* Write a snp block with the given information to the given ostream object.
		* Genotype probabilities must be supplied by the given GenotypeProbabilityGetter
		* objects, which must be callable as
		* - get_AA_probability( index )
		* - get_AB_probability( index )
		* - get_BB_probability( index )
		* where index is the index of the individual in the SNP block.
*/
		template< typename GenotypeProbabilityGetter >
		void write_snp_block(
			std::ostream& aStream,
			uint32_t const flags,
			uint32_t number_of_samples,
			unsigned char max_id_size,
			std::string SNPID,
			std::string RSID,
			unsigned char chromosome,
			uint32_t SNP_position,
			std::string first_allele,
			std::string second_allele,
			GenotypeProbabilityGetter get_AA_probability,
			GenotypeProbabilityGetter get_AB_probability,
			GenotypeProbabilityGetter get_BB_probability
		) ;

		// Implementation

		namespace impl {
			double get_probability_conversion_factor( uint32_t flags ) ;

			template< typename IntegerType >
			void read_length_followed_by_data( std::istream& in_stream, IntegerType* length_ptr, std::string* string_ptr ) {
				IntegerType& length = *length_ptr ;
				read_little_endian_integer( in_stream, length_ptr ) ;
				std::vector< char >buffer ( length ) ;
				in_stream.read( &buffer[0], length ) ;
				string_ptr->assign( buffer.begin(), buffer.end() ) ;
			}

			template< typename IntegerType >
			void write_length_followed_by_data( std::ostream& out_stream, IntegerType length, std::string const data_string ) {
				assert( length <= data_string.size() ) ;
				write_little_endian_integer( out_stream, length ) ;
				out_stream.write( data_string.data(), length ) ;
			}

			template< typename FloatType >
			uint16_t round_to_nearest_integer( FloatType number ) {
				return static_cast< uint16_t > ( number + 0.5 ) ;
			}

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
				
			void read_snp_identifying_data(
				std::istream& aStream,
				uint32_t const flags,
				uint32_t* number_of_samples,
				std::string* SNPID,
				std::string* RSID,
				unsigned char* chromosome,
				uint32_t* SNP_position,
				std::string* first_allele,
				std::string* second_allele
			) ;

			void write_snp_identifying_data(
				std::ostream& aStream,
				uint32_t const flags,
				uint32_t number_of_samples,
				unsigned char max_id_size,
				std::string SNPID,
				std::string RSID,
				unsigned char chromosome,
				uint32_t SNP_position,
				std::string first_allele,
				std::string second_allele
			) ;

			void ignore_snp_probability_data(
				std::istream& aStream,
				uint32_t const flags,
				uint32_t number_of_samples
			) ;

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

			template<
				typename GenotypeProbabilitySetter
			>
			void read_snp_probability_data(
				std::istream& aStream,
				uint32_t const flags,
				uint32_t number_of_samples,
				GenotypeProbabilitySetter set_genotype_probabilities,
				std::vector< char >* buffer1,
				std::vector< char >* buffer2
			) {
				uint32_t const data_size = 6 * number_of_samples ;
				buffer1->resize( data_size ) ;
				if( flags & bgen::e_CompressedSNPBlocks ) {
					uint32_t compressed_data_size ;
					genfile::read_little_endian_integer( aStream, &compressed_data_size ) ;
					buffer2->resize( compressed_data_size ) ;
					aStream.read( &(*buffer2)[0], compressed_data_size ) ;
					zlib_uncompress( *buffer2, buffer1 ) ;
					assert( buffer1->size() == data_size ) ;
				}
				else {
					aStream.read( &(*buffer1)[0], data_size ) ;
				}
				bgen::impl::read_uncompressed_snp_probability_data(
					&(*buffer1)[0],
					&(*buffer1)[0] + data_size,
					impl::get_probability_conversion_factor( flags ),
					number_of_samples,
					set_genotype_probabilities
				) ;
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
			void write_uncompressed_snp_block(
				std::ostream& aStream,
				uint32_t const flags,
				uint32_t const number_of_samples,
				unsigned char max_id_size,
				std::string SNPID,
				std::string RSID,
				unsigned char chromosome,
				uint32_t SNP_position,
				std::string first_allele,
				std::string second_allele,
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

				write_snp_identifying_data( aStream, flags, number_of_samples, max_id_size, SNPID, RSID, chromosome, SNP_position, first_allele, second_allele ) ;
				aStream.write( &uncompressed_buffer[0], uncompressed_data_size ) ;
			}


			template< typename GenotypeProbabilityGetter >
			void write_compressed_snp_block(
				std::ostream& aStream,
				uint32_t const flags,
				uint32_t number_of_samples,
				unsigned char max_id_size,
				std::string SNPID,
				std::string RSID,
				unsigned char chromosome,
				uint32_t SNP_position,
				std::string first_allele,
				std::string second_allele,
				GenotypeProbabilityGetter get_AA_probability,
				GenotypeProbabilityGetter get_AB_probability,
				GenotypeProbabilityGetter get_BB_probability
			) {
				write_snp_identifying_data( aStream, flags, number_of_samples, max_id_size, SNPID, RSID, chromosome, SNP_position, first_allele, second_allele ) ;
				write_compressed_snp_probability_data( aStream, flags, number_of_samples, get_AA_probability, get_AB_probability, get_BB_probability ) ;
			}

			void read_snp_probability_data(
				std::istream& aStream,
				uint32_t const flags,
				uint32_t number_of_samples,
				Ignorer const&
			) ;
		}

		template< typename OffsetType >
		void read_offset( std::istream& iStream, OffsetType* offset ) {
			uint32_t real_offset = 0 ;
			read_little_endian_integer( iStream, &real_offset ) ;
			*offset = real_offset ;
		}

		template< typename OffsetType >
		void write_offset( std::ostream& oStream, OffsetType offset ) {
			uint32_t real_offset = offset ;
			write_little_endian_integer( oStream, real_offset ) ;
		}

		template<
		typename HeaderSizeSetter,
				 typename NumberOfSNPBlocksSetter,
				 typename NumberOfSamplesSetter,
				 typename FlagsSetter,
				 typename FreeDataSetter,
				 typename VersionSetter
		>
		void read_header_block(
			std::istream& aStream,
			HeaderSizeSetter set_header_size,
			NumberOfSNPBlocksSetter set_number_of_snp_blocks,
			NumberOfSamplesSetter set_number_of_samples,
			FreeDataSetter set_free_data,
			FlagsSetter set_flags,
			VersionSetter set_version = VersionSetter()
		) {
			uint32_t
				header_size = 0,
				number_of_snp_blocks = 0,
				number_of_samples = 0,
				flags = 0 ;

			char version[4] ;

			std::size_t fixed_data_size = get_header_block_size( "" ) ;
			std::vector<char> free_data ;

			genfile::read_little_endian_integer( aStream, &header_size ) ;
			assert( header_size >= fixed_data_size ) ;
			genfile::read_little_endian_integer( aStream, &number_of_snp_blocks ) ;
			genfile::read_little_endian_integer( aStream, &number_of_samples ) ;
			aStream.read( &(version[0]), 4 ) ;
			free_data.resize( header_size - fixed_data_size ) ;
			aStream.read( &(free_data[0]), free_data.size() ) ;
			genfile::read_little_endian_integer( aStream, &flags ) ;

			if( aStream ) {
				set_header_size( header_size ) ;
				set_number_of_snp_blocks( number_of_snp_blocks ) ;
				set_number_of_samples( number_of_samples ) ;
				set_free_data( std::string( free_data.begin(), free_data.end() )) ;
				set_flags( flags ) ;
				set_version( std::string( version, version + 4 )) ;
			} else {
				throw BGenError() ;
			}
		}

		template<
			typename NumberOfSamplesSetter,
			typename SampleIdentifierSetter
		>
		void read_sample_block(
			std::istream& aStream,
			NumberOfSamplesSetter set_number_of_samples,
			SampleIdentifierSetter set_identifier
		) {
			uint32_t block_size = 0 ;
			uint32_t number_of_samples ;

			genfile::read_little_endian_integer( aStream, &block_size ) ;
			genfile::read_little_endian_integer( aStream, &number_of_samples ) ;
			set_number_of_samples( number_of_samples ) ;

			uint16_t identifier_size ;
			std::string identifier ;

			for( std::size_t i = 0; i < number_of_samples; ++i ) {
				impl::read_length_followed_by_data( aStream, &identifier_size, &identifier ) ;
				if( aStream ) {
					set_identifier( i, identifier ) ;
				} else {
					throw BGenError() ;
				}
			}
		}

		template< typename GenotypeProbabilityGetter >
		void write_snp_block(
			std::ostream& aStream,
			uint32_t const flags,
			uint32_t const number_of_samples,
			unsigned char max_id_size,
			std::string SNPID,
			std::string RSID,
			unsigned char chromosome,
			uint32_t SNP_position,
			std::string first_allele,
			std::string second_allele,
			GenotypeProbabilityGetter get_AA_probability,
			GenotypeProbabilityGetter get_AB_probability,
			GenotypeProbabilityGetter get_BB_probability
		) {
			if( flags & e_CompressedSNPBlocks ) {
				impl::write_compressed_snp_block(
					aStream,
					flags,
					number_of_samples,
					max_id_size,
					SNPID,
					RSID,
					chromosome,
					SNP_position,
					first_allele,
					second_allele,
					get_AA_probability,
					get_AB_probability,
					get_BB_probability
				) ;
			}
			else {
				impl::write_uncompressed_snp_block(
					aStream,
					flags,
					number_of_samples,
					max_id_size,
					SNPID,
					RSID,
					chromosome,
					SNP_position,
					first_allele,
					second_allele,
					get_AA_probability,
					get_AB_probability,
					get_BB_probability
				) ;
			}
		}
	}
}

#endif
