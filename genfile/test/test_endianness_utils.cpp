#include <iostream>
#include <stdint.h>
#include "test_case.hpp"
#include "genfile/endianness_utils.hpp"

namespace {
	template< typename T, int SizeOf >
	struct LittleEndianTester {
		void operator()( char const* data ) const ;
	} ;

	template< typename T >
	struct LittleEndianTester< T, 1 > {
		void operator()( char const* data, T const expected ) const {
			T t ;
			char const* buffer = data ;
			buffer = genfile::read_little_endian_integer( buffer, buffer + 8, &t ) ;
			TEST_ASSERT( buffer == data + 1 ) ;
			TEST_ASSERT( t == *reinterpret_cast< T const* >( data ) ) ;
			for( int i = 0; i < 1; ++i ) {
				TEST_ASSERT( *( reinterpret_cast< char* >( &t ) + i ) == data[i] ) ;
			}
		}
	} ;

	template< typename T >
	struct LittleEndianTester< T, 2 > {
		void operator()( char const* data, T const expected ) const {
			T t ;
			char const* buffer = data ;
			buffer = genfile::read_little_endian_integer( buffer, buffer + 8, &t ) ;
			TEST_ASSERT( buffer == data + 2 ) ;
			TEST_ASSERT( t == *reinterpret_cast< T const* >( data ) ) ;
			for( int i = 0; i < 2; ++i ) {
				TEST_ASSERT( *( reinterpret_cast< char* >( &t ) + i ) == data[i] ) ;
			}
		}
	} ;

	template< typename T >
	struct LittleEndianTester< T, 4 > {
		void operator()( char const* data, T const expected ) const {
			T t ;
			char const* buffer = data ;
			buffer = genfile::read_little_endian_integer( buffer, buffer + 8, &t ) ;
			TEST_ASSERT( buffer == data + 4 ) ;
			TEST_ASSERT( t == *reinterpret_cast< T const* >( data ) ) ;
			for( int i = 0; i < 4; ++i ) {
				TEST_ASSERT( *( reinterpret_cast< char* >( &t ) + i ) == data[i] ) ;
			}
		}
	} ;

	template< typename T >
	struct LittleEndianTester< T, 8 > {
		void operator()( char const* data, T const expected ) const {
			T t ;
			char const* buffer = data ;
			buffer = genfile::read_little_endian_integer( buffer, buffer + 8, &t ) ;
			TEST_ASSERT( buffer == data + 8 ) ;
			TEST_ASSERT( t == *reinterpret_cast< T const* >( data ) ) ;
			for( int i = 0; i < 8; ++i ) {
				TEST_ASSERT( *( reinterpret_cast< char* >( &t ) + i ) == data[i] ) ;
			}
		}
	} ;
	
	template< typename T >
	LittleEndianTester< T, sizeof(T) > little_endian_tester() {
		return LittleEndianTester< T, sizeof( T ) >() ;
	}

	bool machine_is_little_endian() {
		unsigned int x = 1 ;
		return *reinterpret_cast< char* >( &x ) == 1 ;
	}
	
	template< typename T >
	void read_little_endian_test( char const* data ) {
		// This is if the machine is little-endian
		if( machine_is_little_endian() ) {
			little_endian_tester< T >()( data, *( reinterpret_cast< T const* >( data ))) ;
		}
		else {
			char reversed[8] ;
			std::reverse_copy( &data[0], &data[0] + 8, &reversed[0] ) ;
			little_endian_tester< T >()( data, *( reinterpret_cast< T const* >( reversed ))) ;
		}
	}
	
	template< typename T >
	void write_and_read_little_endian_test( char const* data ) {
		T const* t = reinterpret_cast< T const* >( data ) ;
		// std::cerr << "Writing and reading value: " << *t << ".\n" ;
		char buffer[8] ;
		TEST_ASSERT(
			genfile::write_little_endian_integer( &buffer[0], &buffer[0] + 8, *t ) == ( buffer + sizeof(T) )
		) ;
		T result ;
		TEST_ASSERT( genfile::read_little_endian_integer( &buffer[0], &buffer[0] + 8, &result ) == ( buffer + sizeof( T ) )) ;
		TEST_ASSERT( result == *t ) ;
	}
}

AUTO_TEST_CASE( test_read_little_endian_integer ) {
	std::cerr << "test_read_little_endian_integer()..." ;
	
	char const raw_data[31] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 } ;
	for( int i = 0; i < ( 31 - 8 ); ++i ) {
		char const* data = &raw_data[0] + i ;
		read_little_endian_test< char >( data ) ;
		read_little_endian_test< unsigned char >( data ) ;
		read_little_endian_test< short >( data ) ;
		read_little_endian_test< unsigned short >( data ) ;
		read_little_endian_test< int >( data ) ;
		read_little_endian_test< unsigned int >( data ) ;
		read_little_endian_test< long >( data ) ;
		read_little_endian_test< unsigned long >( data ) ;
		read_little_endian_test< long long >( data ) ;
		read_little_endian_test< unsigned long long >( data ) ;

		read_little_endian_test< int8_t >( data ) ;
		read_little_endian_test< uint8_t >( data ) ;
		read_little_endian_test< int16_t >( data ) ;
		read_little_endian_test< uint16_t >( data ) ;
		read_little_endian_test< int32_t >( data ) ;
		read_little_endian_test< uint32_t >( data ) ;
		read_little_endian_test< int64_t >( data ) ;
		read_little_endian_test< uint64_t >( data ) ;
	}
	std::cerr << "done.\n" ;
}

AUTO_TEST_CASE( test_write_and_read_little_endian_integer ) {
	std::cerr << "test_write_and_read_little_endian_integer()..." ;
	
	char raw_data[265+255] ;
	for( int i = 0; i < 256; ++i ) {
		raw_data[i] = char( i ) ;
	}
	for( int i = 256; i < ( 256 + 255 ); ++i ) {
		raw_data[i] = char( 512 - i - 1 ) ;
	}

	for( int i = 0; i < ( 256 + 255 - 1 ); ++i ) {
		char const* data = &raw_data[0] + i ;
		write_and_read_little_endian_test< char >( data ) ;
		write_and_read_little_endian_test< unsigned char >( data ) ;
		write_and_read_little_endian_test< short >( data ) ;
		write_and_read_little_endian_test< unsigned short >( data ) ;
		write_and_read_little_endian_test< int >( data ) ;
		write_and_read_little_endian_test< unsigned int >( data ) ;
		write_and_read_little_endian_test< long >( data ) ;
		write_and_read_little_endian_test< unsigned long >( data ) ;
		write_and_read_little_endian_test< long long >( data ) ;
		write_and_read_little_endian_test< unsigned long long >( data ) ;

		write_and_read_little_endian_test< int8_t >( data ) ;
		write_and_read_little_endian_test< uint8_t >( data ) ;
		write_and_read_little_endian_test< int16_t >( data ) ;
		write_and_read_little_endian_test< uint16_t >( data ) ;
		write_and_read_little_endian_test< int32_t >( data ) ;
		write_and_read_little_endian_test< uint32_t >( data ) ;
		write_and_read_little_endian_test< int64_t >( data ) ;
		write_and_read_little_endian_test< uint64_t >( data ) ;
	}
	std::cerr << "done.\n" ;
}
