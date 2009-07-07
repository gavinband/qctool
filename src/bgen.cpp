#include <iostream>
#include "endianness_utils.hpp"
#include "bgen.hpp"

namespace gen {
	namespace bgen {
		void write_header_block(
			std::ostream& aStream,
			uint32_t number_of_snp_blocks,
			uint32_t number_of_samples,
			uint32_t snp_block_size,
			std::string const& free_data,
			uint32_t flags
		) {
			unsigned int fixed_data_size
				= sizeof( uint32_t ) + sizeof number_of_snp_blocks + sizeof number_of_samples + sizeof snp_block_size + sizeof flags ; 
			impl::uint32_t header_size = fixed_data_size + free_data.size() ;

			impl::write_little_endian_integer( aStream, header_size ) ;
			impl::write_little_endian_integer( aStream, number_of_snp_blocks ) ;
			impl::write_little_endian_integer( aStream, number_of_samples ) ;
			impl::write_little_endian_integer( aStream, snp_block_size ) ;
			aStream.write( free_data.data(), free_data.size() ) ;
			impl::write_little_endian_integer( aStream, flags ) ;
		}
	}
}