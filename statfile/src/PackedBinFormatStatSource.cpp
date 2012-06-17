
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <string>
#include "statfile/StatSource.hpp"
#include "statfile/IstreamAggregator.hpp"
#include "statfile/PackedBinFormatStatSource.hpp"
#include "statfile/endianness_utils.hpp"

namespace statfile {
	PackedBinFormatStatSource::PackedBinFormatStatSource( std::string const& filename )
	{
		setup( filename ) ;
	}
	
	PackedBinFormatStatSource::PackedBinFormatStatSource( std::auto_ptr< std::istream > stream_ptr )
	{
		setup( stream_ptr ) ;
	}

	std::size_t PackedBinFormatStatSource::number_of_rows() const { return m_number_of_rows ; }
	std::string PackedBinFormatStatSource::get_descriptive_text() const { return m_descriptive_text ; }

	void PackedBinFormatStatSource::reset_to_start() {
		reset_stream_to( m_header_block_length ) ;
		base_t::reset_to_start() ;
		m_next_entry_rows_read = 0 ;
		m_next_entry_column = 0 ;
		get_next_entry_location() ;
	}

	void PackedBinFormatStatSource::get_next_entry_location() {
		uint32_t next_entry_rows_read = 0 ;
		uint32_t next_entry_column = 0 ;
		read_little_endian_integer( stream(), &next_entry_rows_read ) ;
		read_little_endian_integer( stream(), &next_entry_column ) ;
		
		if(
			(next_entry_rows_read < m_next_entry_rows_read)
			|| ((next_entry_rows_read == m_next_entry_rows_read) && ( next_entry_column < m_next_entry_column ) )
		) {
			throw FormatInvalidError() ;
		}

		m_next_entry_rows_read = next_entry_rows_read ;
		m_next_entry_column = next_entry_column ;
	}

	void PackedBinFormatStatSource::read_value( int32_t& value ) {
		if( current_column() == m_next_entry_column && number_of_rows_read() == m_next_entry_rows_read ) {
			read_type( 'i' ) ;
			read_little_endian_integer( stream(), &value ) ;
			get_next_entry_location() ;
		}
		else {
			value = int32_t(0) ;
		}
	}

	void PackedBinFormatStatSource::read_value( uint32_t& value ) {
		if( current_column() == m_next_entry_column && number_of_rows_read() == m_next_entry_rows_read ) {
			read_type( 'u' ) ;
			read_little_endian_integer( stream(), &value ) ;
			get_next_entry_location() ;
		}
		else {
			value = uint32_t(0) ;
		}
	}

	void PackedBinFormatStatSource::read_value( std::string& value ) {
		if( current_column() == m_next_entry_column && number_of_rows_read() == m_next_entry_rows_read ) {
			read_string_value( value ) ;
			get_next_entry_location() ;
		}
		else {
			value = "" ;
		}
	}

	void PackedBinFormatStatSource::read_string_value( std::string& value ) {
		read_type( 's' ) ;
		uint32_t size ;
		read_little_endian_integer( stream(), &size ) ;
		std::vector< char > buffer( size ) ;
		stream().read( &buffer[0], size ) ;
		value.assign( buffer.begin(), buffer.end() ) ;
	}

	void PackedBinFormatStatSource::read_value( double& value ) {
		if( current_column() == m_next_entry_column && number_of_rows_read() == m_next_entry_rows_read ) {
			read_type( 'd' ) ;
			stream().read( reinterpret_cast< char* >( &value ), sizeof( value )) ;
			get_next_entry_location() ;
		} else {
			value = 0.0 ;
		}
	}

	void PackedBinFormatStatSource::ignore_value() {
		if( current_column() == m_next_entry_column && number_of_rows_read() == m_next_entry_rows_read ) {
			ignore_next_value() ;
		}
		else {
			// do nothing, as there's no entry to ignore.
		}
	}

	void PackedBinFormatStatSource::ignore_next_value() {
		char type = read_type() ;
		std::string svalue ;
		switch( type ) {
			case 'i':
			case 'u':
				stream().ignore( 5 ) ;
				break ;
			case 'd':
				stream().ignore( 1 + sizeof( double )) ;
				break ; 
			case 's':
				read_string_value( svalue ) ;
				break ;
			default:
				assert(0) ;
		}
		get_next_entry_location() ;
	}
	
	void PackedBinFormatStatSource::ignore_all() {
		uint32_t number_of_rows_now = number_of_rows_read() ;
		while( m_next_entry_rows_read == number_of_rows_now ) {
			ignore_next_value() ;
		}
	}

	// Read a type specifier, without removing it from the stream.
	char PackedBinFormatStatSource::read_type() {
		char type = stream().get();
		if( stream() && type != 'i' && type != 'u' && type != 'd' && type != 's' ) {
			throw FormatInvalidError() ;
		}
		stream().putback( type ) ;
		return type ;
	}

	void PackedBinFormatStatSource::read_type( char expected ) {
		char type = read_type() ;
		if( stream() ) {
			if( type != expected ) {
				throw FormatInvalidError() ;
			}
			stream().get() ; // remove type specifier from stream.
		}
	}

	void PackedBinFormatStatSource::setup( std::auto_ptr< std::istream > stream_ptr ) {
		set_stream( stream_ptr ) ;
		read_header_block() ;
		m_next_entry_rows_read = 0 ;
		m_next_entry_column = 0 ;
		get_next_entry_location() ;
	}

	void PackedBinFormatStatSource::setup( std::string const& filename ) {
		setup( open_binary_file_for_input( filename )) ;
	}

	void PackedBinFormatStatSource::read_header_block() {
		stream().seekg( 0, std::ios::beg ) ;

		read_type( 'u' ) ;
		read_little_endian_integer( stream(), &m_header_block_length ) ;
		read_string_value( m_version_text ) ;
		assert( m_version_text == "packedbin_v1" ) ;
		read_type( 'u' ) ;
		read_little_endian_integer( stream(), &m_number_of_rows ) ;
		read_type( 'u' ) ;
		uint32_t number_of_columns ;
		read_little_endian_integer( stream(), &number_of_columns ) ;
		read_type( 'u' ) ;
		uint32_t dummy ;
		read_little_endian_integer( stream(), &dummy ) ;
		assert( dummy == number_of_columns ) ;
		read_string_value( m_descriptive_text ) ;
		// now read column headers.
		for( std::size_t i = 0; stream() && i < number_of_columns; ++i ) {
			std::string column_name ;
			read_string_value( column_name ) ;
			if( stream() ) {
				add_column( column_name ) ;
			}
		}
	}
}
