
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_ENDIANNESS_UTILS_HPP
#define GENFILE_ENDIANNESS_UTILS_HPP

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <limits>
#include <cassert>
#include "stdint.h"

namespace genfile {
	// Read an integer stored in little-endian format into an integer stored in memory.
	template< typename IntegerType >
	char const* read_little_endian_integer( char const* buffer, char const* const end, IntegerType* integer_ptr ) {
		assert( end >= buffer + sizeof( IntegerType )) ;
		*integer_ptr = 0 ;
		for( std::size_t byte_i = 0; byte_i < sizeof( IntegerType ); ++byte_i ) {
			(*integer_ptr) |= IntegerType( *reinterpret_cast< unsigned char const* >( buffer++ )) << ( 8 * byte_i ) ;
		}
		return buffer ;
	}

	// Read an integer stored in little-endian format into an integer stored in memory.
	// The stream is assumed to have sizeof( Integertype ) readable bytes.
	template< typename IntegerType >
	void read_little_endian_integer( std::istream& in_stream, IntegerType* integer_ptr ) {
		char buffer[ sizeof( IntegerType ) ] ;
		in_stream.read( buffer, sizeof( IntegerType )) ;
		if( in_stream ) {
			read_little_endian_integer( buffer, buffer + sizeof( IntegerType ), integer_ptr ) ;
		}
	}

	// Read an integer stored in big-endian format into an integer stored in memory.
	template< typename IntegerType >
	char const* read_big_endian_integer( char const* buffer, char const* const end, IntegerType* integer_ptr ) {
		assert( end >= buffer + sizeof( IntegerType )) ;
		*integer_ptr = 0 ;
		for( std::size_t byte_i = 0; byte_i < sizeof( IntegerType ); ++byte_i ) {
			(*integer_ptr) |= IntegerType( *reinterpret_cast< unsigned char const* >( buffer++ )) << ( 8 * ( sizeof( IntegerType ) - byte_i - 1 ) ) ;
		}
		return buffer ;
	}

	// Read an integer stored in big-endian format into an integer stored in memory.
	// The stream is assumed to have sizeof( Integertype ) readable bytes.
	template< typename IntegerType >
	void read_big_endian_integer( std::istream& in_stream, IntegerType* integer_ptr ) {
		char buffer[ sizeof( IntegerType ) ] ;
		in_stream.read( buffer, sizeof( IntegerType )) ;
		if( in_stream ) {
			read_big_endian_integer( buffer, buffer + sizeof( IntegerType ), integer_ptr ) ;
		}
	}

	// Write an integer to the buffer in little-endian format.
	template< typename IntegerType >
	char* write_little_endian_integer( char* buffer, char* const end, IntegerType const integer ) {
		assert( end >= buffer + sizeof( IntegerType )) ;
		for( std::size_t byte_i = 0; byte_i < sizeof( IntegerType ); ++byte_i ) {
			*buffer++ = ( integer >> ( 8 * byte_i ) ) & 0xff ;
		}
		return buffer ;
	}

	// Write an integer to the stream in little-endian format.
	// The stream is assumed to have sizeof( Integertype ) bytes writeable.
	template< typename IntegerType >
	void write_little_endian_integer( std::ostream& out_stream, IntegerType const integer ) {
		char buffer[ sizeof( IntegerType ) ] ;
		write_little_endian_integer( buffer, buffer + sizeof( IntegerType ), integer ) ;
		out_stream.write( buffer, sizeof( IntegerType )) ;
	}

	// Write an integer to the buffer in big-endian format.
	template< typename IntegerType >
	char* write_big_endian_integer( char* buffer, char* const end, IntegerType const integer ) {
		assert( end >= buffer + sizeof( IntegerType )) ;
		for( std::size_t byte_i = 0; byte_i < sizeof( IntegerType ); ++byte_i ) {
			*buffer++ = ( integer >> ( 8 * ( sizeof( IntegerType ) - byte_i - 1 ) ) ) & 0xff ;
		}
		return buffer ;
	}

	// Write an integer to the stream in big-endian format.
	// The stream is assumed to have sizeof( Integertype ) bytes writeable.
	template< typename IntegerType >
	void write_big_endian_integer( std::ostream& out_stream, IntegerType const integer ) {
		char buffer[ sizeof( IntegerType ) ] ;
		write_big_endian_integer( buffer, buffer + sizeof( IntegerType ), integer ) ;
		out_stream.write( buffer, sizeof( IntegerType )) ;
	}

	template< typename IntegerType >
	void read_length_followed_by_data( std::istream& in_stream, IntegerType* length_ptr, std::string* string_ptr ) {
		IntegerType& length = *length_ptr ;
		read_little_endian_integer( in_stream, length_ptr ) ;
		std::vector< char > buffer( length ) ;
		in_stream.read( &buffer[0], length ) ;
		string_ptr->assign( &buffer[0], &buffer[0] + length ) ;
	}

	template< typename IntegerType >
	void write_length_followed_by_data( std::ostream& out_stream, IntegerType length, std::string const data_string ) {
		assert( length <= data_string.size() ) ;
		write_little_endian_integer( out_stream, length ) ;
		out_stream.write( data_string.data(), length ) ;
	}

	namespace impl {
		template< typename IntegerType, bool is_signed >
		struct small_integer_writer
		{
		} ;

		template< typename IntegerType >
		struct small_integer_writer< IntegerType, true >
		{
			char* write( char* buffer, char* const end, IntegerType const integer ) const {
				assert( end >= buffer + 1 ) ;
				if( integer == 0 ) {
					*buffer++ = char( 0 | 0x80 )  ;
				}
				else if( integer == int8_t( integer )) {
					*buffer++ = char( sizeof( uint8_t ) | 0x80 ) ;
					buffer = write_little_endian_integer( buffer, end, int8_t( integer )) ;
				}
				else if( integer == int16_t( integer )) {
					*buffer++ = char( sizeof( uint16_t ) | 0x80 ) ;
					buffer = write_little_endian_integer( buffer, end, int16_t( integer )) ;
				}
				else if( integer == int32_t( integer )) {
					*buffer++ = char( sizeof( uint32_t ) | 0x80 ) ;
					buffer = write_little_endian_integer( buffer, end, int32_t( integer )) ;
				}
				else if( integer == int64_t( integer )) {
					*buffer++ = char( sizeof( uint64_t ) | 0x80 ) ;
					buffer = write_little_endian_integer( buffer, end, int64_t( integer )) ;
				}
				else {
					assert(0) ;
				}
				return buffer ;
			}

			std::size_t size( IntegerType const integer ) const {
				if( integer == 0 ) {
					return 1 ;
				}
				else if( integer == int8_t( integer )) {
					return 2 ; 
				}
				else if( integer == int16_t( integer )) {
					return 2 ; 
				}
				else if( integer == int32_t( integer )) {
					return 4 ; 
				}
				else if( integer == int64_t( integer )) {
					return 8 ; 
				}
				else {
					assert(0) ;
				}
			}

			char const* read( char const* buffer, char const* const end, IntegerType* integer ) const {
				assert( end >= buffer + 1 ) ;
				char num_bytes = *buffer++ ;
				assert( num_bytes & 0x80 ) ; // only allow reading from a signed integer type.
				num_bytes = num_bytes & 0x7F ;
				if( num_bytes == 0 ) {
					*integer = 0 ;
				}
				else if( num_bytes == 1 ) {
					int8_t value ;
					buffer = read_little_endian_integer( buffer, end, &value ) ;
					*integer = value ;
					assert( *integer == value ) ;
				}
				else if( num_bytes == 2 ) {
					int16_t value ;
					buffer = read_little_endian_integer( buffer, end, &value ) ;
					*integer = value ;
					assert( *integer == value ) ;
				}
				else if( num_bytes == 4 ) {
					int32_t value ;
					buffer = read_little_endian_integer( buffer, end, &value ) ;
					*integer = value ;
					assert( *integer == value ) ;
				}
				else if( num_bytes == 8 ) {
					int64_t value ;
					buffer = read_little_endian_integer( buffer, end, &value ) ;
					*integer = value ;
					assert( *integer == value ) ;
				}
				else {
					assert(0) ;
				}
				return buffer ;
			}
		} ;

		template< typename IntegerType >
		struct small_integer_writer< IntegerType, false >
		{
			char* write( char* buffer, char* const end, IntegerType const integer ) const {
				assert( end >= buffer + 1 ) ;
				if( integer == 0 ) {
					*buffer++ = char( 0 )  ;
				}
				else if( integer == uint8_t( integer )) {
					*buffer++ = char( sizeof( uint8_t ) ) ;
					buffer = write_little_endian_integer( buffer, end, uint8_t( integer )) ;
				}
				else if( integer == uint16_t( integer )) {
					*buffer++ = char( sizeof( uint16_t )) ;
					buffer = write_little_endian_integer( buffer, end, uint16_t( integer )) ;
				}
				else if( integer == uint32_t( integer )) {
					*buffer++ = char( sizeof( uint32_t ) ) ;
					buffer = write_little_endian_integer( buffer, end, uint32_t( integer )) ;
				}
				else if( integer == uint64_t( integer )) {
					*buffer++ = char( sizeof( uint64_t ) ) ;
					buffer = write_little_endian_integer( buffer, end, uint64_t( integer )) ;
				}
				else {
					assert(0) ;
				}
				return buffer ;
			}
			

			std::size_t size( IntegerType const integer ) const {
				if( integer == 0 ) {
					return 1 ;
				}
				if( integer < std::numeric_limits< uint8_t >::max() ) {
					return 2 ;
				}
				else if( integer < std::numeric_limits< uint16_t >::max() ) {
					return 3 ;
				}
				else if( integer < std::numeric_limits< uint32_t >::max() ) {
					return 4 ;
				}
				else if( integer < std::numeric_limits< uint64_t >::max() ) {
					return 8 ;
				}
				else {
					assert(0) ;
				}
			}
			
			char const* read( char const* buffer, char const* const end, IntegerType* integer ) const {
				assert( end >= buffer + 1 ) ;
				char num_bytes = *buffer++ ;
				assert(( num_bytes & 0x80 ) == 0 ) ; // only allow reading into an unsigned integer type.
				if( num_bytes == 0 ) {
					*integer = 0 ;
				}
				else if( num_bytes == 1 ) {
					uint8_t value ;
					buffer = read_little_endian_integer( buffer, end, &value ) ;
					*integer = value ;
					assert( *integer == value ) ;
				}
				else if( num_bytes == 2 ) {
					uint16_t value ;
					buffer = read_little_endian_integer( buffer, end, &value ) ;
					*integer = value ;
					assert( *integer == value ) ;
				}
				else if( num_bytes == 4 ) {
					uint32_t value ;
					buffer = read_little_endian_integer( buffer, end, &value ) ;
					*integer = value ;
					assert( *integer == value ) ;
				}
				else if( num_bytes == 8 ) {
					uint64_t value ;
					buffer = read_little_endian_integer( buffer, end, &value ) ;
					*integer = value ;
					assert( *integer == value ) ;
				}
				else {
					assert(0) ;
				}
				return buffer ;
			}
		} ;
	}

	// Write an integer to the buffer, fitting it into the smallest number of bytes
	// that will hold it.  A single character is written first to record the type.
	template< typename IntegerType >
	char* write_small_integer( char* buffer, char* const end, IntegerType const integer ) {
		return impl::small_integer_writer< IntegerType, std::numeric_limits< IntegerType >::is_signed >().write( buffer, end, integer ) ;
	}

	// Get the number of bytes write_small_integer will use to store an integer.
	template< typename IntegerType >
	std::size_t get_small_integer_size( IntegerType const integer ) {
		return impl::small_integer_writer< IntegerType, std::numeric_limits< IntegerType >::is_signed >().size( integer ) ;
	}
	
	// Read an integer from the buffer, in the representation of write_small_integer().
	template< typename IntegerType >
	char const* read_small_integer( char const* buffer, char const* const end, IntegerType* integer ) {
		return impl::small_integer_writer< IntegerType, std::numeric_limits< IntegerType >::is_signed >().read( buffer, end, integer ) ;
	}
}

#endif

