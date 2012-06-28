
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
		namespace impl {
			typedef ::uint32_t uint32_t ;
			typedef ::uint16_t uint16_t ;
		}

		typedef impl::uint32_t uint32_t ;
		typedef impl::uint16_t uint16_t ;
		enum FlagsType { e_NoFlags = 0, e_CompressedSNPBlocks = 0x1, e_MultiCharacterAlleles = 0x2, e_LongIds = 0x4 } ;

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

		template<
			typename IntegerSetter,
			typename StringSetter,
			typename AlleleSetter,
			typename ChromosomeSetter,
			typename SNPPositionSetter,
			typename GenotypeProbabilitySetter
		>
		void read_compressed_snp_block(
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
		) ;

		// Implementation

		namespace impl {
			double const PROBABILITY_CONVERSION_FACTOR = 10000.0 ;

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

			template< typename IntegerType >
			void read_little_endian_integer( std::istream& in_stream, IntegerType* integer_ptr ) {
				genfile::read_little_endian_integer( in_stream, integer_ptr ) ;
			}

			template< typename IntegerType >
			void write_little_endian_integer( std::ostream& out_stream, IntegerType integer ) {
				genfile::write_little_endian_integer( out_stream, integer ) ;
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

			template<
				typename GenotypeProbabilitySetter
			>
			void read_snp_probability_data(
				std::istream& aStream,
				uint32_t const flags,
				uint32_t number_of_samples,
				GenotypeProbabilitySetter set_genotype_probabilities
			) {
				for ( uint32_t i = 0 ; i < number_of_samples ; ++i ) {
					uint16_t AA = 0, AB = 0, BB = 0 ;
					read_little_endian_integer( aStream, &AA ) ;
					read_little_endian_integer( aStream, &AB ) ;
					read_little_endian_integer( aStream, &BB ) ;

					double const factor = ( bool( flags & e_LongIds ) ) ? 32768.0 : PROBABILITY_CONVERSION_FACTOR ;
					set_genotype_probabilities(
						i,
						convert_from_integer_representation( AA, factor ),
						convert_from_integer_representation( AB, factor ),
						convert_from_integer_representation( BB, factor )
					) ;
				}
			}

			void read_snp_probability_data(
				std::istream& aStream,
				uint32_t const flags,
				uint32_t number_of_samples,
				Ignorer const&
			) ;

			template<
				typename GenotypeProbabilitySetter
			>
			void read_compressed_snp_probability_data(
				std::istream& aStream,
				uint32_t const flags,
				uint32_t number_of_samples,
				GenotypeProbabilitySetter set_genotype_probabilities
			) {
				// read the size of the compressed data
				uint32_t compressed_data_size = 0 ;
				impl::read_little_endian_integer( aStream, &compressed_data_size ) ;
				uLongf uncompressed_data_size = (6 * number_of_samples) ;
				// Construct buffers for the uncompressed and compressed data
				std::vector< Bytef > compressed_data_buffer( compressed_data_size ) ;
				std::vector< Bytef > uncompressed_data_buffer( uncompressed_data_size ) ;
				// Read the data and uncompress
				aStream.read( reinterpret_cast< char*> (&compressed_data_buffer[0]), compressed_data_size ) ;
				int result = uncompress( &uncompressed_data_buffer[0], &uncompressed_data_size, &compressed_data_buffer[0], compressed_data_size ) ;
				assert( result == Z_OK ) ;
				// Load the uncompressed data into an istringstream
				std::istringstream sStream ;
				sStream.rdbuf()->pubsetbuf( reinterpret_cast< char* >( &uncompressed_data_buffer[0] ), uncompressed_data_size ) ;

				if ( aStream ) {
					// If all ok thus far, we'll take the plunge, calling the callbacks to (presumably)
					// set up the user's data structure.  This way we avoid allocating memory
					// for the genotype probabilities here.	 Note that if an error occurs while reading
					// the probabilities, the user's data may therefore be left in a state which does
					// not correspond to any actual valid snp block.
					read_snp_probability_data( sStream, flags, number_of_samples, set_genotype_probabilities ) ;
				}
			}

			void read_compressed_snp_probability_data(
				std::istream& aStream,
				uint32_t const flags,
				uint32_t number_of_samples,
				Ignorer const&
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

			template< typename GenotypeProbabilityGetter >
			void write_snp_probability_data(
				std::ostream& aStream,
				uint32_t const flags,
				uint32_t number_of_samples,
				GenotypeProbabilityGetter get_AA_probability,
				GenotypeProbabilityGetter get_AB_probability,
				GenotypeProbabilityGetter get_BB_probability
			) {
				double const factor = ( bool( flags & e_LongIds ) ? 32768.0 : impl::PROBABILITY_CONVERSION_FACTOR ) ;
				for ( impl::uint32_t i = 0 ; i < number_of_samples ; ++i ) {
					impl::uint16_t
						AA = convert_to_integer_representation( get_AA_probability( i ), factor ),
						AB = convert_to_integer_representation( get_AB_probability( i ), factor ),
						BB = convert_to_integer_representation( get_BB_probability( i ), factor ) ;

					write_little_endian_integer( aStream, AA ) ;
					write_little_endian_integer( aStream, AB ) ;
					write_little_endian_integer( aStream, BB ) ;
				}
			}
		}

		template< typename OffsetType >
		void read_offset( std::istream& iStream, OffsetType* offset ) {
			impl::uint32_t real_offset = 0 ;
			read_little_endian_integer( iStream, &real_offset ) ;
			*offset = real_offset ;
		}

		template< typename OffsetType >
		void write_offset( std::ostream& oStream, OffsetType offset ) {
			impl::uint32_t real_offset = offset ;
			write_little_endian_integer( oStream, real_offset ) ;
		}

		template<
		typename HeaderSizeSetter,
		typename NumberOfSNPBlocksSetter,
		typename NumberOfSamplesSetter,
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
		) {
			impl::uint32_t
			header_size = 0,
			number_of_snp_blocks = 0,
			number_of_samples = 0,
			reserved = 0,
			flags = 0 ;

			std::size_t fixed_data_size
			= get_header_block_size( "" ) ;

			std::vector<char> free_data ;

			impl::read_little_endian_integer( aStream, &header_size ) ;
			assert( header_size >= fixed_data_size ) ;
			impl::read_little_endian_integer( aStream, &number_of_snp_blocks ) ;
			impl::read_little_endian_integer( aStream, &number_of_samples ) ;
			impl::read_little_endian_integer( aStream, &reserved ) ;
			free_data.resize( header_size - fixed_data_size ) ;
			aStream.read( &(free_data[0]), header_size - fixed_data_size ) ;
			impl::read_little_endian_integer( aStream, &flags ) ;

			if ( aStream ) {
				set_header_size( header_size ) ;
				set_number_of_snp_blocks( number_of_snp_blocks ) ;
				set_number_of_samples( number_of_samples ) ;
				set_free_data( std::string( free_data.begin(), free_data.end() )) ;
				set_flags( flags ) ;
			}
		}

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
		) {
			// We read the following data in the following order.
			// Initialisers are provided so that, in the case where the stream becomes bad (e.g. end of file)
			// the assertions below are not triggered.
			impl::uint32_t number_of_samples = 0;
			std::string SNPID ;
			std::string RSID ;
			unsigned char chromosome ;
			impl::uint32_t SNP_position = 0;
			std::string first_allele, second_allele ;

			impl::read_snp_identifying_data( aStream, flags, &number_of_samples, &SNPID, &RSID, &chromosome, &SNP_position, &first_allele, &second_allele ) ;

			if ( aStream ) {
				// If all ok thus far, we'll take the plunge, calling the callbacks to (presumably)
				// set up the user's data structure.  This way we avoid allocating memory
				// for the genotype probabilities here.	 Note that if an error occurs while reading
				// the probabilities, the user's data may therefore be left in a state which does
				// not correspond to any actual valid snp block.
				set_number_of_samples( number_of_samples ) ;
				set_SNPID( SNPID ) ;
				set_RSID( RSID ) ;
				set_chromosome( chromosome ) ;
				set_SNP_position( SNP_position ) ;
				set_allele1( first_allele ) ;
				set_allele2( second_allele ) ;

				impl::read_snp_probability_data( aStream, flags, number_of_samples, set_genotype_probabilities ) ;
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
			impl::write_snp_identifying_data( aStream, flags, number_of_samples, max_id_size, SNPID, RSID, chromosome, SNP_position, first_allele, second_allele ) ;
			impl::write_snp_probability_data( aStream, flags, number_of_samples, get_AA_probability, get_AB_probability, get_BB_probability ) ;
		}


		template<
			typename IntegerSetter,
			typename StringSetter,
			typename AlleleSetter,
			typename ChromosomeSetter,
			typename SNPPositionSetter,
			typename GenotypeProbabilitySetter
		>
		void read_compressed_snp_block(
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
		) {
#if !HAVE_ZLIB
			assert(0) ; // zlib is required for compression support.
#else
			// We read the following data in the following order.
			// Initialisers are provided so that, in the case where the stream becomes bad (e.g. end of file)
			// the assertions below are not triggered.
			impl::uint32_t number_of_samples = 0;
			std::string SNPID ;
			std::string RSID ;
			unsigned char chromosome ;
			impl::uint32_t SNP_position = 0;
			std::string first_allele, second_allele ;

			impl::read_snp_identifying_data( aStream, flags, &number_of_samples, &SNPID, &RSID, &chromosome, &SNP_position, &first_allele, &second_allele ) ;

			if ( aStream ) {
				set_number_of_samples( number_of_samples ) ;
				set_SNPID( SNPID ) ;
				set_RSID( RSID ) ;
				set_chromosome( chromosome ) ;
				set_SNP_position( SNP_position ) ;
				set_allele1( first_allele ) ;
				set_allele2( second_allele ) ;

				impl::read_compressed_snp_probability_data( aStream, flags, number_of_samples, set_genotype_probabilities ) ;
			}
#endif
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
#if !HAVE_ZLIB
			assert(0) ; // zlib is required for compression support.
#else
			impl::write_snp_identifying_data( aStream, flags, number_of_samples, max_id_size, SNPID, RSID, chromosome, SNP_position, first_allele, second_allele ) ;

			// Construct a buffer into which we will compress.
			uLongf uncompressed_data_size = (6 * number_of_samples) ;
			uLongf buffer_size = 12 + (1.1 * uncompressed_data_size) ;		// calculated according to zlib manual.
			std::vector< Bytef > compression_buffer( buffer_size ) ;
			std::vector< Bytef > uncompressed_buffer( uncompressed_data_size ) ;
			// get probability data, uncompressed.
			std::ostringstream oStream ;
			oStream.rdbuf()->pubsetbuf( reinterpret_cast< char* >( &uncompressed_buffer[0] ), uncompressed_data_size ) ;
			impl::write_snp_probability_data( oStream, flags, number_of_samples, get_AA_probability, get_AB_probability, get_BB_probability ) ;
			// compress it
			int result = compress( &compression_buffer[0], &buffer_size, &uncompressed_buffer[0], uncompressed_data_size ) ;
			assert( result == Z_OK ) ;
			// and write it (buffer_size is now the compressed length of the data).
			uint32_t the_buffer_size = buffer_size ;
			impl::write_little_endian_integer( aStream, the_buffer_size ) ;
			aStream.write( reinterpret_cast< char const* >( &compression_buffer[0] ), buffer_size ) ;
#endif
		}
	}
}

#endif
