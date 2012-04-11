
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <string>
#include "statfile/StatSink.hpp"
#include "statfile/OstreamAggregator.hpp"
#include "statfile/endianness_utils.hpp"
#include "statfile/BinFormatStatSink.hpp"

namespace statfile {
	// Construct from a filename
	BinFormatStatSink::BinFormatStatSink( std::string const& filename, uint32_t number_of_id_columns )
		: m_number_of_id_columns( number_of_id_columns )
	{
		setup( filename ) ;
	}

	// Construct from a stream (must be opened in binary mode).
	BinFormatStatSink::BinFormatStatSink( std::auto_ptr< std::ostream > stream_ptr , uint32_t number_of_id_columns )
		: m_number_of_id_columns( number_of_id_columns )
	{
		setup( stream_ptr ) ;
	}

	// Finalise written file
	BinFormatStatSink::~BinFormatStatSink() {
		write_header_block() ;
	}

	void BinFormatStatSink::set_descriptive_text( std::string const& descriptive_text ) {
		// ensure we haven't written any data...
		assert( number_of_rows_written() == 0 && current_column() == 0 ) ;
		m_descriptive_text = descriptive_text ;
		write_header_block() ;
	}

	void BinFormatStatSink::write_value( int32_t const& value ) {
		assert( current_column() < m_number_of_id_columns ) ;
		stream().put( 'i' ) ;
		write_little_endian_integer( stream(), value ) ;		
	}

	void BinFormatStatSink::write_value( uint32_t const& value ) {
		assert( current_column() < m_number_of_id_columns ) ;
		stream().put( 'u' ) ;
		write_little_endian_integer( stream(), value ) ;		
	}

	void BinFormatStatSink::write_value( std::string const& value ) {
		assert( current_column() < m_number_of_id_columns ) ;
		stream().put( 's' ) ;
		uint32_t size = value.size() ;
		write_little_endian_integer( stream(), size ) ;
		stream().write( value.data(), size ) ;
	}

	void BinFormatStatSink::write_value( double const& value ) {
		if( current_column() < m_number_of_id_columns ) {
			stream().put( 'd' ) ;
		}
		stream().write( reinterpret_cast< char const* >( &value ), sizeof( value )) ;
	}

	void BinFormatStatSink::add_column_impl( std::string const& name ) {
		base_t::add_column_impl( name ) ;
		write_header_block() ;
	}

	void BinFormatStatSink::setup( std::auto_ptr< std::ostream > stream_ptr ) {
		set_stream( stream_ptr ) ;
		m_version_text = "statfilebin_v1" ;
		m_descriptive_text = "" ;
		write_header_block() ;
	}

	void BinFormatStatSink::setup( std::string const& filename ) {
		setup( open_binary_file_for_output( filename )) ;
	}

	// This function must match write_header_block below.
	uint32_t BinFormatStatSink::get_header_block_length() const { 
		uint32_t result = 
			5									// header length
			+ ( m_version_text.size() + 5 )		// version text
			+ 5									// number of rows
			+ 5									// number of columns
			+ 5									// number of id columns.
			+ ( m_descriptive_text.size() + 5 ) // descriptive text
		;

		for( std::size_t i = 0; i < number_of_columns(); ++i ) {
			result += ( column_name( i ).size() + 5 ) ;
		}

		return result ;
	}
	
	// Write the header block -- its length must equal get_header_length().
	void BinFormatStatSink::write_header_block() {
		// position after offset bytes.
		stream().seekp( 0, std::ios::beg ) ;
		uint32_t const number_of_rows = number_of_rows_written() ;
		uint32_t const no_of_columns = number_of_columns() ;

		stream().put( 'u' ) ;
		write_little_endian_integer( stream(), get_header_block_length() ) ;		
		stream().put( 's' ) ;
		write_little_endian_integer( stream(), uint32_t( m_version_text.size() )) ;		
		stream().write( m_version_text.data(), m_version_text.size() ) ;
		stream().put( 'u' ) ;
		write_little_endian_integer( stream(), number_of_rows ) ;		
		stream().put( 'u' ) ;
		write_little_endian_integer( stream(), no_of_columns ) ;		
		stream().put( 'u' ) ;
		write_little_endian_integer( stream(), m_number_of_id_columns ) ;		
		stream().put( 's' ) ;
		write_little_endian_integer( stream(), uint32_t( m_descriptive_text.size() )) ;		
		stream().write( m_descriptive_text.data(), m_descriptive_text.size() ) ;

		// Now write column names.
		for( std::size_t i = 0; i < no_of_columns; ++i ) {
			stream().put( 's' ) ;
			write_little_endian_integer( stream(), uint32_t( column_name( i ).size() )) ;		
			stream().write( column_name( i ).data(), column_name( i ).size() ) ;
		}
		stream().flush() ;
	}
}
