#include <iostream>
#include <functional>
#include <cassert>
#include "endianness_utils.hpp"
#include "genbin.hpp"

namespace genfilebin {
	namespace impl {
		uint16_t round_to_nearest_integer( double number ) {
			return static_cast< uint16_t > ( number + 0.5 ) ;
		}
	}

	void read_header_block( std::istream& aStream, HeaderBlock* header_block ) {
		read_header_block(
			aStream,
			std::bind1st( std::mem_fun( &HeaderBlock::set_number_of_individuals ), header_block ),
			std::bind1st( std::mem_fun( &HeaderBlock::set_ID_field_storage ), header_block ),
			std::bind1st( std::mem_fun( &HeaderBlock::set_free_data ), header_block ),
			std::bind1st( std::mem_fun( &HeaderBlock::set_number_of_blocks ), header_block )
		) ;
	}

	void write_header_block(
		std::ostream& aStream,
		HeaderBlock const& header_block
	) {
		return write_header_block(
			aStream,
			header_block.number_of_individuals(),
			header_block.ID_field_storage(),
			header_block.free_data(),
			header_block.number_of_blocks()
		) ;
	}
	
	void write_header_block(
		std::ostream& aStream,
		impl::uint32_t number_of_individuals,
		unsigned char ID_field_storage,
		std::string free_data,
		impl::uint32_t number_of_blocks
	) {
		impl::uint32_t flags = 0 ;
		impl::uint32_t header_length = 13 + free_data.size() ;

		assert( number_of_blocks == 1 ) ;  // As per the spec.
		impl::write_little_endian_integer( aStream, header_length ) ;
		impl::write_little_endian_integer( aStream, number_of_individuals ) ;
		impl::write_little_endian_integer( aStream, flags ) ;
		impl::write_little_endian_integer( aStream, ID_field_storage ) ;
		aStream.write( free_data.data(), free_data.size() ) ;
		impl::write_little_endian_integer( aStream, number_of_blocks ) ;
	}
}