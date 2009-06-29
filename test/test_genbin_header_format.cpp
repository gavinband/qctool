#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include "../config.hpp"
#if HAVE_BOOST_UNIT_TEST
	#define BOOST_AUTO_TEST_MAIN
	#include "boost/test/auto_unit_test.hpp"
	#define AUTO_TEST_CASE( param ) BOOST_AUTO_TEST_CASE(param)
	#define TEST_ASSERT( param ) BOOST_ASSERT( param )
#else
	#define AUTO_TEST_CASE( param ) void param()
	#define TEST_ASSERT( param ) assert( param )
#endif
#include "GenRow.hpp"
#include "genbin.hpp"
#include "stdint.h"

namespace data {
	std::string construct_header_block(
		uint32_t number_of_samples,
		unsigned char ID_field_storage,
		std::string free_data,
		uint32_t number_of_blocks
	) {
		std::ostringstream oStream ;
		uint32_t header_length = free_data.size() + 13 ;
		write_little_endian_integer( oStream, header_length ) ;
		write_little_endian_integer( oStream, number_of_samples ) ;
		write_little_endian_integer( oStream, ID_field_storage ) ;
		oStream.write( free_data.data(), free_data.size() ) ;
		write_little_endian_integer( oStream, number_of_blocks ) ;

		return oStream.str() ;
	}

	std::string to_hex( std::string const& str ) {
		std::ostringstream o ;
		for( std::size_t i = 0; i < str.size(); ++i ) {
			if( i % 4 == 0 )
				o << "|" ;
			o << std::hex << std::setw(2) << std::setfill('0') << static_cast<int> (str[i]) ;
		}
		return o.str() ;
	}
}

template< typename T >
struct Setter
{
	Setter( T& field ): m_field( field ) {} ;
	void operator()( T const& value ) { m_field = value ; }
private:
	T& m_field ;
} ;

template< typename T >
Setter< T > make_setter( T& field ) { return Setter<T>( field ) ; }

void do_header_block_read_test(
	uint32_t number_of_individuals,
	unsigned char ID_field_storage,
	std::string free_data,
	uint32_t number_of_blocks
) {
	std::istringstream inStream ;
	inStream.str(
		data::construct_header_block(
			number_of_individuals,
			ID_field_storage,
			free_data,
			number_of_blocks
		)
	) ;

	uint32_t number_of_individuals2 ;
	unsigned char ID_field_storage2 ;
	std::string free_data2 ;
	uint32_t number_of_blocks2 ;

	genbin::read_header_block(
		inStream,
		make_setter( number_of_individuals2 ),
		make_setter( ID_field_storage2 ),
		make_setter( free_data2 ),
		make_setter( number_of_blocks2 )
	) ;
	
	TEST_ASSERT( inStream ) ;
	TEST_ASSERT( number_of_individuals == number_of_individuals2 ) ;	
	TEST_ASSERT( ID_field_storage == ID_field_storage2 ) ;	
	TEST_ASSERT( free_data == free_data2 ) ;
	TEST_ASSERT( number_of_blocks == number_of_blocks2 ) ;	
}

void do_header_block_write_test( 
	uint32_t number_of_individuals,
	unsigned char ID_field_storage,
	std::string free_data,
	uint32_t number_of_blocks
) {
	std::ostringstream outStream ;
	genbin::write_header_block( 
		outStream,
		number_of_individuals,
		ID_field_storage,
		free_data,
		number_of_blocks
	) ;

	std::string expected = data::construct_header_block(
		number_of_individuals,
		ID_field_storage,
		free_data,
		number_of_blocks
	) ;
	
	std::cout << "\"" << data::to_hex( outStream.str() ) << "\"\n" ;
	std::cout << "\"" << data::to_hex( expected ) << "\"\n" ;

	TEST_ASSERT( outStream.str() == expected ) ;
}

AUTO_TEST_CASE( test_header_block_input ) {
	std::cout << "test_header_block_input\n" ;
	do_header_block_read_test( 100, 45, "Hi, this is some data", 1 ) ;
	do_header_block_read_test( 4294967295u, 255, "Ochen ochen, ochen ura. coheh   heheh  he  agg agg 767  $$%%$   **  ", 4294967295u ) ;
	do_header_block_read_test( 0, 0, "", 1 ) ;
	do_header_block_read_test( 0, 0, "", 0 ) ;
}

AUTO_TEST_CASE( test_header_block_output ) {
	std::cout << "test_header_block_output\n" ;
	do_header_block_write_test( 6, 6, "This is some free text.", 1 ) ;
	do_header_block_write_test( 4294967295u, 25, "Some more free text here, you see &&&& $$$ $$$", 1 ) ;
	do_header_block_write_test( 0, 0, "", 1 ) ;
}

#ifndef HAVE_BOOST_UNIT_TEST

int main( int argc, char** argv ) {
	test_header_block_input() ;
	test_header_block_output() ;
}

#endif

