#ifndef GENBIN_REFERENCE_IMPLEMENTATION_HPP
#define GENBIN_REFERENCE_IMPLEMENTATION_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <stdint.h>
#include "endianness_utils.hpp"

/*
* This file contains a reference implementation of the GENBIN file format
* specification described at:
* http://www.well.ox.ac.uk/~gav/genbin_format.html
*
* LICENSE: you are free to use this for personal and/or academic purposes only.
*
* To use this file you will also need "endianness_utils.hpp".
*
* This implementation provides the following functions.
*
* - genbin::read_header_block( stream, data... )
* - genbin::write_header_block( stream, data... )
* - genbin::read_snp_block( stream, data... )
* - genbin::write_snp_block( stream, data... )
*
*/

namespace genfilebin {
	namespace impl {
		typedef ::uint32_t uint32_t ;
		typedef ::uint16_t uint16_t ;
		template< typename T > void ignore( T const& ) {} ;
	}

	struct HeaderBlock
	{
		public:
			HeaderBlock()
				: m_number_of_individuals( 0 ),
				  m_ID_field_storage( 50 ),
				  m_free_data( "" ),
				  m_number_of_blocks(1)
			{}
			
		impl::uint32_t number_of_individuals() const { return m_number_of_individuals ; }
		unsigned char ID_field_storage() const { return m_ID_field_storage ; }
		std::string const& free_data() const { return m_free_data ; }
		impl::uint32_t number_of_blocks() const { return m_number_of_blocks ; }
		impl::uint32_t header_length() const { return 13 + m_free_data.size() ; }

		HeaderBlock& set_number_of_individuals( impl::uint32_t number_of_individuals ) { m_number_of_individuals = number_of_individuals ; return *this ; }
		HeaderBlock& set_ID_field_storage( unsigned char ID_field_storage ) { m_ID_field_storage = ID_field_storage ; return *this ; }
		HeaderBlock& set_free_data( std::string free_data ) { m_free_data = free_data ; return *this ; }
		HeaderBlock& set_number_of_blocks( impl::uint32_t number_of_blocks ) { m_number_of_blocks = number_of_blocks ; return *this ; }

		private:
			impl::uint32_t m_number_of_individuals ;
			unsigned char m_ID_field_storage ;
			std::string m_free_data ;
			impl::uint32_t m_number_of_blocks ;
	} ;

	/*
	* Function: read_offset
	* Read the offset from the start of the stream.
	*/
	template< typename OffsetType >
	void read_offset( std::istream& iStream, impl::uint32_t* offset ) ;

	/*
	* Function: read_offset
	* Read the offset value to the stream.
	*/
	template< typename OffsetType >
	void write_offset( std::ostream& oStream, impl::uint32_t offset ) ;

	/*
	* Function: read_header_block()
	* Read a header block from the given istream object, returning the data in a HeaderBlock struct.
	*/
	void read_header_block( std::istream& aStream, HeaderBlock* header_block ) ;

	/*
	* Function: write_header_block()
	* Write a header block with the given information to the given ostream object.
	* Returns: the number of bytes written.
	*/
	void write_header_block(
		std::ostream& aStream,
		HeaderBlock const& header_block
	) ;

	/*
	* Function: read_header_block()
	* Read a header block from the given istream object, returning the data using
	* the setter objects (or function pointers) passed as arguments.
	*/
	template<
		typename NumberOfIndividualsSetter,
		typename IDFieldStorageSetter,
		typename FreeDataSetter,
		typename NumberOfBlocksSetter >
	void read_header_block(
		std::istream& aStream,
		NumberOfIndividualsSetter set_number_of_individuals,
		IDFieldStorageSetter set_ID_field_storage,		
		FreeDataSetter set_free_data,
		NumberOfBlocksSetter set_number_of_blocks = &impl::ignore< impl::uint32_t >
	) ;

	/*
	* Function: write_header_block()
	* Write a header block with the given information to the given ostream object.
	* Note: currently, the value of number_of_blocks must be 1.
	* Returns: the number of bytes written.
	*/
	void write_header_block(
		std::ostream& aStream,
		impl::uint32_t number_of_individuals,
		unsigned char ID_field_storage,
		std::string free_data,
		impl::uint32_t number_of_blocks
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
		typename AllelesSetter,
		typename SNPPositionSetter,
		typename GenotypeProbabilitySetter
	>
	void read_snp_block(
		std::istream& aStream,
		IntegerSetter set_number_of_samples,
		StringSetter set_SNPID,
		StringSetter set_RSID,
		SNPPositionSetter set_SNP_position,
		AllelesSetter set_alleles,
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
		impl::uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		impl::uint32_t SNP_position,
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
		
		uint16_t round_to_nearest_integer( double number ) ;
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
		typename NumberOfIndividualsSetter,
		typename IDFieldStorageSetter,
		typename FreeDataSetter,
		typename NumberOfBlocksSetter
	>
	void read_header_block(
		std::istream& aStream,
		NumberOfIndividualsSetter set_number_of_individuals,
		IDFieldStorageSetter set_ID_field_storage,		
		FreeDataSetter set_free_data,
		NumberOfBlocksSetter set_number_of_blocks = &(impl::ignore< impl::uint32_t >)
	) {
		impl::uint32_t
			header_length,
			number_of_blocks,
			number_of_individuals,		
			flags ;
		unsigned char S ;
		std::vector<char> free_data ;

		impl::read_little_endian_integer( aStream, &header_length ) ;
		assert( header_length >= 13 ) ;
		impl::read_little_endian_integer( aStream, &number_of_individuals ) ;
		impl::read_little_endian_integer( aStream, &flags ) ;
		assert( flags == 0 ) ;
		S = aStream.get() ;
		free_data.resize( header_length - 13 ) ;
		aStream.read( &free_data[0], header_length - 13 ) ;
		impl::read_little_endian_integer( aStream, &number_of_blocks ) ;
		 // assert( number_of_blocks == 1 ) ;
		
		if( aStream ) {
			set_number_of_individuals( number_of_individuals ) ;
			set_ID_field_storage( S ) ;
			set_free_data( std::string( free_data.begin(), free_data.end() )) ;
			set_number_of_blocks( number_of_blocks ) ;
		}
	}

	template<
		typename IntegerSetter,
		typename StringSetter,
		typename AllelesSetter,
		typename SNPPositionSetter,
		typename GenotypeProbabilitySetter
	>
	void read_snp_block(
		std::istream& aStream,
		IntegerSetter set_number_of_samples,
		StringSetter set_SNPID,
		StringSetter set_RSID,
		SNPPositionSetter set_SNP_position,
		AllelesSetter set_alleles,
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
			set_alleles( first_allele, second_allele ) ;

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
		impl::uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		impl::uint32_t SNP_position,
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

#endif
