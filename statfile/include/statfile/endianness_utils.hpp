#ifndef STATFILE_ENDIAN_HPP
#define STATFILE_ENDIAN_HPP

#include <iostream>
#include <algorithm>
#include <cassert>
#include "stdint.h"

namespace statfile {
	
	namespace impl {
		uint32_t const i = 0x04030201 ;
		char const* const ptr = reinterpret_cast< char const * const >( &i ) ;
		bool const compile_machine_is_little_endian = ( *(ptr+0) == 1 ) && ( *(ptr+1) == 2) && ( *(ptr+2) == 3 ) && ( *(ptr+3) == 4 );
		bool const compile_machine_is_big_endian = ( *(ptr+0) == 4 ) && ( *(ptr+1) == 3) && ( *(ptr+2) == 2 ) && ( *(ptr+3) == 1 );
	}

	// Read an integer stored in little-endian format into an integer stored
	// in memory, taking account of the machine's endianness.
	// The stream is assumed to have as many bytes readable as the integer's in-memory representation.
	template< typename IntegerType >
	void read_little_endian_integer( std::istream& in_stream, IntegerType* integer_ptr ) {
		char* ptr = reinterpret_cast< char* >( integer_ptr ) ;
		in_stream.read( ptr, sizeof( IntegerType )) ;
		if( impl::compile_machine_is_little_endian ) {
			// do nothing
		}
		else if( impl::compile_machine_is_big_endian ) {
			std::reverse( ptr, ptr + sizeof( IntegerType )) ;
		}
		else {
			assert( 0 ) ; // only little and big endian supported.
		}
	}

	// Write an integer to the stream in little-endian format.
	// The stream is assumed to have as many bytes writeable as the integer's in-memory representation.
	template< typename IntegerType >
	void write_little_endian_integer( std::ostream& out_stream, IntegerType const integer ) {
		char const* ptr = reinterpret_cast< char const* >( &integer ) ;
		if( impl::compile_machine_is_little_endian ) {
			out_stream.write( ptr, sizeof( IntegerType )) ;
		}
		else if( impl::compile_machine_is_big_endian ) {
			IntegerType buffer ;
			char* buffer_ptr = reinterpret_cast< char* >( &buffer ) ;
			std::reverse_copy( ptr, ptr + sizeof( IntegerType ), buffer_ptr ) ;
			out_stream.write( ptr, sizeof( IntegerType )) ;
		}
		else {
			assert( 0 ) ; // Only little and big endian supported.
		}
	}
}

#endif