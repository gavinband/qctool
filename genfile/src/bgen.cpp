
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include "genfile/Chromosome.hpp"
#include "genfile/endianness_utils.hpp"
#include "genfile/bgen.hpp"

namespace genfile {
	namespace bgen {
		namespace impl {
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
			) {

				if( aStream ) {
					impl::read_little_endian_integer( aStream, number_of_samples ) ;
				}
				
				if( flags & e_LongIds ) {
					uint16_t SNPID_size = 0;
					uint16_t RSID_size = 0;
					uint32_t allele1_size = 0;
					uint32_t allele2_size = 0;
					uint16_t chromosome_size = 0 ;
					std::string chromosome_string ;
					impl::read_length_followed_by_data( aStream, &SNPID_size, SNPID ) ;
					impl::read_length_followed_by_data( aStream, &RSID_size, RSID ) ;
					impl::read_length_followed_by_data( aStream, &chromosome_size, &chromosome_string ) ;
					impl::read_little_endian_integer( aStream, SNP_position ) ;
					impl::read_length_followed_by_data( aStream, &allele1_size, first_allele ) ;
					impl::read_length_followed_by_data( aStream, &allele2_size, second_allele ) ;
					*chromosome = ChromosomeEnum( Chromosome( chromosome_string ) ) ;
				}
				else {
					unsigned char max_id_size = 0;
					unsigned char SNPID_size = 0;
					unsigned char RSID_size = 0;
					
					if( aStream ) {
						impl::read_little_endian_integer( aStream, &max_id_size ) ;
					}
					if( aStream ) {
						impl::read_length_followed_by_data( aStream, &SNPID_size, SNPID ) ;
						assert( SNPID_size <= max_id_size ) ;
						aStream.ignore( max_id_size - SNPID_size ) ;
					}
					if( aStream ) {
						impl::read_length_followed_by_data( aStream, &RSID_size, RSID ) ;
						assert( RSID_size <= max_id_size ) ;
						aStream.ignore( max_id_size - RSID_size ) ;
					}
					if( aStream ) {
						impl::read_little_endian_integer( aStream, chromosome ) ;
						impl::read_little_endian_integer( aStream, SNP_position ) ;

						if( flags & e_MultiCharacterAlleles ) {
							unsigned char allele1_size = 0 ;
							unsigned char allele2_size = 0 ;
							impl::read_length_followed_by_data( aStream, &allele1_size, first_allele ) ;
							assert( allele1_size <= max_id_size ) ;
							impl::read_length_followed_by_data( aStream, &allele2_size, second_allele ) ;
							assert( allele2_size <= max_id_size ) ;
						} else {
							first_allele->resize( 1 ) ;
							(*first_allele)[0] = aStream.get() ;
							second_allele->resize( 1 ) ;
							(*second_allele)[0] = aStream.get() ;
						}
					}
				}
			}
			
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
			) {
				write_little_endian_integer( aStream, number_of_samples ) ;

				if( flags & e_LongIds ) {
					std::size_t const max_allele_length = std::numeric_limits< uint32_t >::max() ;
					std::size_t const max_id_length = std::numeric_limits< uint16_t >::max() ;
					assert( SNPID.size() <= static_cast< std::size_t >( max_id_length )) ;
					assert( RSID.size() <= static_cast< std::size_t >( max_id_length )) ;
					if( first_allele.size() > static_cast< std::size_t >( max_allele_length ) ) {
						std::cerr << "Warning: at SNP " << SNPID << " " << RSID << " pos=" << SNP_position << ", truncating first allele of size " << first_allele.size() << ".\n" ;
						first_allele.resize( max_allele_length - 3 ) ;
						first_allele += "..." ;
					}
					if( second_allele.size() > static_cast< std::size_t >( max_allele_length ) ) {
						std::cerr << "Warning: at SNP " << SNPID << " " << RSID << " pos=" << SNP_position << ", truncating second allele of size " << second_allele.size() << ".\n" ;
						second_allele.resize( max_allele_length - 3 ) ;
						second_allele += "..." ;
					}
					write_length_followed_by_data( aStream, uint16_t( SNPID.size() ), SNPID.data() ) ;
					write_length_followed_by_data( aStream, uint16_t( RSID.size() ), RSID.data() ) ;
					std::string const chromosome_string = Chromosome( chromosome ) ;
					write_length_followed_by_data( aStream, uint16_t( chromosome_string.size() ), chromosome_string ) ;
					write_little_endian_integer( aStream, SNP_position ) ;
					write_length_followed_by_data( aStream, uint32_t( first_allele.size() ), first_allele.data() ) ;
					write_length_followed_by_data( aStream, uint32_t( second_allele.size() ), second_allele.data() ) ;
				}
				else {
					assert( SNPID.size() <= static_cast< std::size_t >( max_id_size )) ;
					assert( RSID.size() <= static_cast< std::size_t >( max_id_size )) ;
					unsigned char SNPID_size = SNPID.size() ;
					unsigned char RSID_size = RSID.size() ;
					SNPID.resize( max_id_size, ' ' ) ;
					RSID.resize( max_id_size, ' ' ) ;

					write_little_endian_integer( aStream, max_id_size ) ;
					write_length_followed_by_data( aStream, SNPID_size, SNPID.data() ) ;
					aStream.write( SNPID.data() + SNPID_size, max_id_size - SNPID_size ) ;
					write_length_followed_by_data( aStream, RSID_size, RSID.data() ) ;
					aStream.write( RSID.data() + RSID_size, max_id_size - RSID_size ) ;
					write_little_endian_integer( aStream, chromosome ) ;
					write_little_endian_integer( aStream, SNP_position ) ;

					if( flags & e_MultiCharacterAlleles ) {
						assert(
							first_allele.size() <= static_cast< std::size_t >( max_id_size )
							&& second_allele.size() <= static_cast< std::size_t >( max_id_size )
						) ;
						unsigned char allele1_size = first_allele.size() ;
						unsigned char allele2_size = second_allele.size() ;
						write_length_followed_by_data( aStream, allele1_size, first_allele.data() ) ;
						write_length_followed_by_data( aStream, allele2_size, second_allele.data() ) ;
					} else {
						assert( first_allele.size() == 1 && second_allele.size() == 1 ) ;
						aStream.put( first_allele[0] ) ;
						aStream.put( second_allele[0] ) ;
					}
				}
			}

			void read_snp_probability_data(
				std::istream& aStream,
				uint32_t const flags,
				uint32_t number_of_samples,
				Ignorer const&
			) {
				aStream.ignore( 6 * number_of_samples ) ;
			}

			void read_compressed_snp_probability_data(
				std::istream& aStream,
				uint32_t const flags,
				uint32_t number_of_samples,
				Ignorer const&
			) {
				// read the size of the compressed data
				uint32_t compressed_data_size ;
				impl::read_little_endian_integer( aStream, &compressed_data_size ) ;
				aStream.ignore( compressed_data_size ) ;
			}
		}
		
		/* Function: get_header_blocK-size()
		* Return the size in bytes of the header block with the given free data
		*/
		std::size_t get_header_block_size(
			std::string const& free_data
		) {
			std::size_t fixed_data_size = 5 * sizeof( uint32_t ) ;
			return fixed_data_size + free_data.size() ;
		}
		
		void write_header_block(
			std::ostream& aStream,
			uint32_t number_of_snp_blocks,
			uint32_t number_of_samples,
			std::string const& free_data,
			uint32_t flags
		) {
			uint32_t reserved = 0u ;
			impl::uint32_t header_size = get_header_block_size( free_data ) ;

			impl::write_little_endian_integer( aStream, header_size ) ;
			impl::write_little_endian_integer( aStream, number_of_snp_blocks ) ;
			impl::write_little_endian_integer( aStream, number_of_samples ) ;
			impl::write_little_endian_integer( aStream, reserved ) ;
			aStream.write( free_data.data(), free_data.size() ) ;
			impl::write_little_endian_integer( aStream, flags ) ;
		}
	}
}
