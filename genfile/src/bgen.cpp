#include <iostream>
#include "endianness_utils.hpp"
#include "bgen.hpp"

namespace genfile {
	namespace bgen {
		namespace impl {
			void read_snp_identifying_data(
                std::istream& aStream,
                uint32_t* number_of_samples,
                std::string* SNPID,
                std::string* RSID,
                uint32_t* SNP_position,
                char* first_allele,
                char* second_allele
            ) {
                unsigned char max_id_size = 0;
                unsigned char SNPID_size = 0;
                unsigned char RSID_size = 0;

                impl::read_little_endian_integer( aStream, number_of_samples ) ;
                impl::read_little_endian_integer( aStream, &max_id_size ) ;
                assert(( impl::MAX_ID_SIZE == 255 ) || ( max_id_size <= impl::MAX_ID_SIZE )) ; // 1st condition avoids warning for MAX_ID_SIZE=255.

                impl::read_length_followed_by_data( aStream, &SNPID_size, SNPID ) ;
                assert( SNPID_size <= max_id_size ) ;
                aStream.ignore( max_id_size - SNPID_size ) ;

                impl::read_length_followed_by_data( aStream, &RSID_size, RSID ) ;
                assert( RSID_size <= max_id_size ) ;
                aStream.ignore( max_id_size - RSID_size ) ;

                impl::read_little_endian_integer( aStream, SNP_position ) ;

                *first_allele = aStream.get() ;
                *second_allele = aStream.get() ;
            }
            
			void write_snp_identifying_data(
                std::ostream& aStream,
                uint32_t number_of_samples,
                unsigned char max_id_size,
                std::string SNPID,
                std::string RSID,
                uint32_t SNP_position,
                char first_allele,
                char second_allele
            ) {
                assert( max_id_size <= impl::MAX_ID_SIZE ) ;
                if ( max_id_size == 0 ) {
                    max_id_size = impl::MAX_ID_SIZE ;
                }

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
                write_little_endian_integer( aStream, SNP_position ) ;
                aStream.put( first_allele ) ;
                aStream.put( second_allele ) ;
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
			uint32_t snp_block_size,
			std::string const& free_data,
			uint32_t flags
		) {
			impl::uint32_t header_size = get_header_block_size( free_data ) ;

			impl::write_little_endian_integer( aStream, header_size ) ;
			impl::write_little_endian_integer( aStream, number_of_snp_blocks ) ;
			impl::write_little_endian_integer( aStream, number_of_samples ) ;
			impl::write_little_endian_integer( aStream, snp_block_size ) ;
			aStream.write( free_data.data(), free_data.size() ) ;
			impl::write_little_endian_integer( aStream, flags ) ;
		}
	}
}