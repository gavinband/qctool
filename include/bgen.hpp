#ifndef GENBIN_REFERENCE_IMPLEMENTATION_HPP
#define GENBIN_REFERENCE_IMPLEMENTATION_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <stdint.h>
#include "endianness_utils.hpp"
#include "snp_data_utils.hpp"

/*
* This file contains a reference implementation of the GENBIN file format
* specification described at:
* http://www.well.ox.ac.uk/~gav/bgen_format.html
*
* LICENSE: you are free to use this for personal and/or academic purposes only.
*
* To use this file you will also need "endianness_utils.hpp".
*
* This implementation provides the following functions.
*
* - bgen::read_offset( stream, data... )
* - bgen::write_offset( stream, data... )
* - bgen::read_snp_block( stream, data... )
* - bgen::write_snp_block( stream, data... )
*/

namespace gen {
	namespace bgen {
		namespace impl {
			typedef ::uint32_t uint32_t ;
			typedef ::uint16_t uint16_t ;
		}
	
		typedef impl::uint32_t uint32_t ;
		typedef impl::uint16_t uint16_t ;

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
			SNPBlockSizeSetter set_snp_block_size,
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
			uint32_t snp_block_size,
			std::string const& free_data,
			uint32_t flags = 0
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
		* - set_alleles( char, char )
		* - set_genotype_probabilities( double, double, double )
		*/
		template<
			typename IntegerSetter,
			typename StringSetter,
			typename AlleleSetter,
			typename SNPPositionSetter,
			typename GenotypeProbabilitySetter
		>
		void read_snp_block(
			std::istream& aStream,
			IntegerSetter set_number_of_samples,
			StringSetter set_SNPID,
			StringSetter set_RSID,
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
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			uint32_t SNP_position,
			char first_allele,
			char second_allele,
			GenotypeProbabilityGetter get_AA_probability,
			GenotypeProbabilityGetter get_AB_probability,
			GenotypeProbabilityGetter get_BB_probability,
			unsigned char max_id_size = 0
		) ;

	// Implementation

		namespace impl {
			unsigned char const MAX_ID_SIZE = 254 ;
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
				::read_little_endian_integer( in_stream, integer_ptr ) ;
			}

			template< typename IntegerType >
			void write_little_endian_integer( std::ostream& out_stream, IntegerType integer ) {
				::write_little_endian_integer( out_stream, integer ) ;
			}

			template< typename FloatType >
			uint16_t round_to_nearest_integer( FloatType number ) {
				return static_cast< uint16_t > ( number + 0.5 ) ;
			}
		}

		template< typename OffsetType >
		void read_offset( std::istream& iStream, OffsetType* offset ) {
			impl::uint32_t real_offset ;
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
			typename SNPBlockSizeSetter,
			typename FlagsSetter,
			typename FreeDataSetter
		>
		void read_header_block(
			std::istream& aStream,
			HeaderSizeSetter set_header_size,
			NumberOfSNPBlocksSetter set_number_of_snp_blocks,
			NumberOfSamplesSetter set_number_of_samples,
			SNPBlockSizeSetter set_snp_block_size,
			FreeDataSetter set_free_data,
			FlagsSetter set_flags
		) {
			impl::uint32_t
				header_size,
				number_of_snp_blocks,
				number_of_samples,		
				snp_block_size,
				flags ;

			unsigned int fixed_data_size
				= sizeof header_size + sizeof number_of_snp_blocks + sizeof number_of_samples + sizeof snp_block_size + sizeof flags ; 

			std::vector<char> free_data ;

			impl::read_little_endian_integer( aStream, &header_size ) ;
			assert( header_size >= fixed_data_size ) ;
			impl::read_little_endian_integer( aStream, &number_of_snp_blocks ) ;
			impl::read_little_endian_integer( aStream, &number_of_samples ) ;
			impl::read_little_endian_integer( aStream, &snp_block_size ) ;
			free_data.resize( header_size - fixed_data_size ) ;
			aStream.read( &free_data[0], header_size - fixed_data_size ) ;
			impl::read_little_endian_integer( aStream, &flags ) ;
			assert( flags == 0 ) ;

			if( aStream ) {
				set_header_size( header_size ) ;
				set_number_of_snp_blocks( number_of_snp_blocks ) ;
				set_number_of_samples( number_of_samples ) ;
				set_snp_block_size( snp_block_size ) ;
				set_free_data( std::string( free_data.begin(), free_data.end() )) ;
			}
		}

		template<
			typename IntegerSetter,
			typename StringSetter,
			typename AlleleSetter,
			typename SNPPositionSetter,
			typename GenotypeProbabilitySetter
		>
		void read_snp_block(
			std::istream& aStream,
			IntegerSetter set_number_of_samples,
			StringSetter set_SNPID,
			StringSetter set_RSID,
			SNPPositionSetter set_SNP_position,
			AlleleSetter set_allele1,
			AlleleSetter set_allele2,
			GenotypeProbabilitySetter set_genotype_probabilities
		) {
			// We read the following data in the following order.
			// Initialisers are provided so that, in the case where the stream becomes bad (e.g. end of file)
			// the assertions below are not triggered.
			impl::uint32_t number_of_samples = 0;
			unsigned char max_id_size = 0;
			unsigned char SNPID_size = 0;
			std::string SNPID ;
			unsigned char RSID_size = 0;	
			std::string RSID ;
			impl::uint32_t SNP_position = 0;
			char first_allele, second_allele ;

			impl::read_little_endian_integer( aStream, &number_of_samples ) ;
			impl::read_little_endian_integer( aStream, &max_id_size ) ;
			assert(( impl::MAX_ID_SIZE == 255 ) || ( max_id_size <= impl::MAX_ID_SIZE )) ; // 1st condition avoids warning for MAX_ID_SIZE=255.

			impl::read_length_followed_by_data( aStream, &SNPID_size, &SNPID ) ;
			assert( SNPID_size <= max_id_size ) ;
			aStream.ignore( max_id_size - SNPID_size ) ;

			impl::read_length_followed_by_data( aStream, &RSID_size, &RSID ) ;
			assert( RSID_size <= max_id_size ) ;
			aStream.ignore( max_id_size - RSID_size ) ;

			impl::read_little_endian_integer( aStream, &SNP_position ) ;

			first_allele = aStream.get() ;
			second_allele = aStream.get() ;

			if( aStream ) {
			// If all ok thus far, we'll take the plunge, calling the callbacks to (presumably)
			// set up the user's data structure.  This way we avoid allocating memory
			// for the genotype probabilities here.  Note that if an error occurs while reading
			// the probabilities, the user's data may therefore be left in a state which does
			// not correspond to any actual valid snp block. 
				set_number_of_samples( number_of_samples ) ;
				set_SNPID( SNPID ) ;
				set_RSID( RSID ) ;
				set_SNP_position( SNP_position ) ;
				set_allele1( first_allele ) ;
				set_allele2( second_allele ) ;

				for( impl::uint32_t i = 0 ; i < number_of_samples ; ++i ) {
					impl::uint16_t AA, AB, BB ;
					impl::read_little_endian_integer( aStream, &AA ) ;
					impl::read_little_endian_integer( aStream, &AB ) ;
					impl::read_little_endian_integer( aStream, &BB ) ;

					set_genotype_probabilities(
						i,
						static_cast< float >( AA ) / impl::PROBABILITY_CONVERSION_FACTOR,
						static_cast< float >( AB ) / impl::PROBABILITY_CONVERSION_FACTOR,
						static_cast< float >( BB ) / impl::PROBABILITY_CONVERSION_FACTOR
					) ;
				}
			}
		}


		template< typename GenotypeProbabilityGetter >
		void write_snp_block(
			std::ostream& aStream,
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			uint32_t SNP_position,
			char first_allele,
			char second_allele,
			GenotypeProbabilityGetter get_AA_probability,
			GenotypeProbabilityGetter get_AB_probability,
			GenotypeProbabilityGetter get_BB_probability,
			unsigned char max_id_size = 0
		) {
			assert( max_id_size <= impl::MAX_ID_SIZE ) ;
			if( max_id_size == 0 ) {
				max_id_size = impl::MAX_ID_SIZE ;
			}

			assert( SNPID.size() <= static_cast< std::size_t >( max_id_size )) ;
			assert( RSID.size() <= static_cast< std::size_t >( max_id_size )) ;
			unsigned char SNPID_size = SNPID.size() ;
			unsigned char RSID_size = RSID.size() ;
			SNPID.resize( max_id_size, ' ' ) ;
			RSID.resize( max_id_size, ' ' ) ;
		
			impl::write_little_endian_integer( aStream, number_of_samples ) ;
			impl::write_little_endian_integer( aStream, max_id_size ) ;
			impl::write_length_followed_by_data( aStream, SNPID_size, SNPID.data() ) ;
			aStream.write( SNPID.data() + SNPID_size, max_id_size - SNPID_size ) ;
			impl::write_length_followed_by_data( aStream, RSID_size, RSID.data() ) ;
			aStream.write( RSID.data() + RSID_size, max_id_size - RSID_size ) ;
			impl::write_little_endian_integer( aStream, SNP_position ) ;
			aStream.put( first_allele ) ;
			aStream.put( second_allele ) ;

			for( impl::uint32_t i = 0 ; i < number_of_samples ; ++i ) {
				impl::uint16_t
					AA = impl::round_to_nearest_integer( get_AA_probability( i ) * impl::PROBABILITY_CONVERSION_FACTOR ),
					AB = impl::round_to_nearest_integer( get_AB_probability( i ) * impl::PROBABILITY_CONVERSION_FACTOR ),
					BB = impl::round_to_nearest_integer( get_BB_probability( i ) * impl::PROBABILITY_CONVERSION_FACTOR ) ;

				impl::write_little_endian_integer( aStream, AA ) ;
				impl::write_little_endian_integer( aStream, AB ) ;
				impl::write_little_endian_integer( aStream, BB ) ;
			}
		}
	}
}

#endif
