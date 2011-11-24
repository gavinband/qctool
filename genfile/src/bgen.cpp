#include <iostream>
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
                unsigned char max_id_size = 0;
                unsigned char SNPID_size = 0;
                unsigned char RSID_size = 0;

				if( aStream ) {
                	impl::read_little_endian_integer( aStream, number_of_samples ) ;
				}
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
                assert( SNPID.size() <= static_cast< std::size_t >( max_id_size )) ;
                assert( RSID.size() <= static_cast< std::size_t >( max_id_size )) ;
                unsigned char SNPID_size = SNPID.size() ;
                unsigned char RSID_size = RSID.size() ;
                SNPID.resize( max_id_size, ' ' ) ;
                RSID.resize( max_id_size, ' ' ) ;

                write_little_endian_integer( aStream, number_of_samples ) ;
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