#include <vector>
#include <string>
#include "statfile/StatSource.hpp"
#include "statfile/IstreamAggregator.hpp"
#include "statfile/BinFormatStatSource.hpp"
#include "statfile/endianness_utils.hpp"

namespace statfile {
	BinFormatStatSource::BinFormatStatSource( std::string const& filename )
	: m_reading_header( false )
	{
		setup( filename ) ;
	}
	
	BinFormatStatSource::BinFormatStatSource( std::auto_ptr< std::istream > stream_ptr )
		: m_reading_header( false )
	{
		setup( stream_ptr ) ;
	}

	void BinFormatStatSource::setup( std::string const& filename ) {
		setup( open_binary_file_for_input( filename )) ;
	}

	void BinFormatStatSource::setup( std::auto_ptr< std::istream > stream_ptr ) {
		set_stream( stream_ptr ) ;
		read_header_block() ;
	}

	void BinFormatStatSource::read_header_block() {
		m_reading_header = true ;
		stream().seekg( 0, std::ios::beg ) ;

		read_value( m_header_block_length ) ;
		read_value( m_version_text ) ;
		assert( m_version_text == "statfilebin_v1" ) ;
		read_value( m_number_of_rows ) ;
		uint32_t number_of_columns ;
		read_value( number_of_columns ) ;
		read_value( m_number_of_id_columns ) ;
		assert( m_number_of_id_columns <= number_of_columns ) ;
		read_value( m_descriptive_text ) ;
		// now read column headers.
		for( std::size_t i = 0; stream() && i < number_of_columns; ++i ) {
			std::string column_name ;
			read_value( column_name ) ;
			if( stream() ) {
				add_column( column_name ) ;
			}
		}
		
		m_reading_header = false ;
	}

	std::size_t BinFormatStatSource::number_of_rows() const { return m_number_of_rows ; }
	std::string const& BinFormatStatSource::get_descriptive_text() const { return m_descriptive_text ; }

	void BinFormatStatSource::reset_to_start() {
		reset_stream_to( m_header_block_length ) ;
		base_t::reset_to_start() ;
	}

	void BinFormatStatSource::read_value( int32_t& value ) {
		assert( m_reading_header || current_column() < m_number_of_id_columns ) ;
		read_type( 'i' ) ;
		read_little_endian_integer( stream(), &value ) ;
	}

	void BinFormatStatSource::read_value( uint32_t& value ) {
		assert( m_reading_header || current_column() < m_number_of_id_columns ) ;
		read_type( 'u' ) ;
		read_little_endian_integer( stream(), &value ) ;
	}

	void BinFormatStatSource::read_value( std::string& value ) {
		assert( m_reading_header || current_column() < m_number_of_id_columns ) ;
		read_type( 's' ) ;
		uint32_t size ;
		read_little_endian_integer( stream(), &size ) ;
		std::vector< char > buffer( size ) ;
		stream().read( &buffer[0], size ) ;
		value.assign( buffer.begin(), buffer.end() ) ;
	}

	void BinFormatStatSource::read_value( double& value ) {
		if( m_reading_header || current_column() < m_number_of_id_columns ) {
			read_type( 'd' ) ;
		}
		stream().read( reinterpret_cast< char* >( &value ), sizeof( value )) ;
	}

	void BinFormatStatSource::ignore_value() {
		if( current_column() < m_number_of_id_columns ) {
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
					read_value( svalue ) ;
					break ;
				default:
					assert(0) ;
			}
		}
		else {
			// value is a double
			stream().ignore(sizeof( double )) ;
		}
	}

	void BinFormatStatSource::ignore_all() {
		std::size_t c = current_column() ;
		for( ; c < m_number_of_id_columns; ++c ) {
			ignore_value() ;		
		}
		stream().ignore((number_of_columns() - c) * sizeof( double )) ;
	}

	// Read a type specifier, without removing it from the stream.
	char BinFormatStatSource::read_type() {
		char type = stream().get();
		if( stream() && type != 'i' && type != 'u' && type != 'd' && type != 's' ) {
			throw FormatInvalidError() ;
		}
		stream().putback( type ) ;
		return type ;
	}

	void BinFormatStatSource::read_type( char expected ) {
		char type = read_type() ;
		if( stream() ) {
			if( type != expected ) {
				throw FormatInvalidError() ;
			}
			stream().get() ; // remove type specifier from stream.
		}
	}

}
