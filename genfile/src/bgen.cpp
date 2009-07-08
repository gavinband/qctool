#include <iostream>
#include "endianness_utils.hpp"
#include "bgen.hpp"

namespace genfile {
	namespace bgen {
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