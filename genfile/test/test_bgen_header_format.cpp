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
#include "genfile/bgen.hpp"
#include "stdint.h"

namespace data {
	std::string construct_header_block(
		uint32_t number_of_snp_blocks,
		uint32_t number_of_samples,
		std::string free_data,
		uint32_t flags
	) {
		uint32_t reserved = 0;
		std::ostringstream oStream ;
		uint32_t header_length = free_data.size() + 20 ;
		write_little_endian_integer( oStream, header_length ) ;
		write_little_endian_integer( oStream, number_of_snp_blocks ) ;
		write_little_endian_integer( oStream, number_of_samples ) ;
		write_little_endian_integer( oStream, reserved ) ;
		oStream.write( free_data.data(), free_data.size() ) ;
		write_little_endian_integer( oStream, flags ) ;

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
	uint32_t number_of_snp_blocks,
	uint32_t number_of_samples,
	std::string free_data,
	uint32_t flags
) {
	std::istringstream inStream ;
	inStream.str(
		data::construct_header_block(
			number_of_snp_blocks,
			number_of_samples,
			free_data,
			flags
		)
	) ;

	uint32_t header_size, number_of_snp_blocks2, number_of_samples2, flags2 ;
	std::string free_data2 ;

	genfile::bgen::read_header_block(
		inStream,
		make_setter( header_size ),
		make_setter( number_of_snp_blocks2 ),
		make_setter( number_of_samples2 ),
		make_setter( free_data2 ),
		make_setter( flags2 )
	) ;
	
	TEST_ASSERT( inStream ) ;
	TEST_ASSERT( header_size == 20 + free_data2.size()) ;
	TEST_ASSERT( number_of_snp_blocks2 == number_of_snp_blocks ) ;	
	TEST_ASSERT( number_of_samples2 == number_of_samples ) ;	
	TEST_ASSERT( free_data == free_data2 ) ;
	TEST_ASSERT( flags2 == flags ) ;	
}

void do_header_block_write_test( 
	uint32_t number_of_snp_blocks,
	uint32_t number_of_samples,
	std::string free_data,
	uint32_t flags
) {
	std::ostringstream outStream ;
	genfile::bgen::write_header_block( 
		outStream,
		number_of_snp_blocks,
		number_of_samples,
		free_data,
		flags
	) ;

	std::string expected = data::construct_header_block(
		number_of_snp_blocks,
		number_of_samples,
		free_data,
		flags
	) ;
	
	std::cout << "\"" << data::to_hex( outStream.str() ) << "\"\n" ;
	std::cout << "\"" << data::to_hex( expected ) << "\"\n" ;

	TEST_ASSERT( outStream.str() == expected ) ;
}

AUTO_TEST_CASE( test_header_block_input ) {
	std::cout << "test_header_block_input\n" ;
	do_header_block_read_test( 100, 45, "Hi, this is some data", 1 ) ;
	do_header_block_read_test( 4294967295u, 4294967295u, "Ochen ochen, ochen ura. coheh   heheh  he  agg agg 767  $$%%$   **  ", 4294967295u ) ;
	do_header_block_read_test( 0, 0, "", 1 ) ;
	do_header_block_read_test( 0, 0, "", 0 ) ;
}

AUTO_TEST_CASE( test_header_block_output ) {
	std::cout << "test_header_block_output\n" ;
	do_header_block_write_test( 6, 6, "This is some free text.", 1 ) ;
	do_header_block_write_test( 4294967295u, 4294967295u, "Some more free text here, you see &&&& $$$ $$$", 4294967295u ) ;
	do_header_block_write_test( 0, 0, "", 1 ) ;
}

#ifndef HAVE_BOOST_UNIT_TEST

int main( int argc, char** argv ) {
	test_header_block_input() ;
	test_header_block_output() ;
}

#endif

