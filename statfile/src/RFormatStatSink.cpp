#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <limits>
#include "statfile/RFormatStatSink.hpp"

namespace statfile {
	RFormatStatSink::RFormatStatSink( std::auto_ptr< std::ostream > stream_ptr )
	: 	m_comment_character( '#' ),
		have_written_header( false ),
		m_precision(6)
	{
		setup( stream_ptr ) ;
	}

	RFormatStatSink::RFormatStatSink( std::string const& filename )
	: 	m_comment_character( '#' ),
		have_written_header( false ),
		m_precision(6)
	{
		setup( filename ) ;
	}

	void RFormatStatSink::set_descriptive_text( std::string const& text ) {
		m_descriptive_text = "#" ;
		for( std::string::const_iterator i = text.begin(); i != text.end(); ++i ) {
			m_descriptive_text += *i ;
			if( *i == '\n' ) {
				m_descriptive_text.append( "#" ) ;
			}
		}
	}

	void RFormatStatSink::setup( std::string const& filename ) {
		setup( open_text_file_for_output( filename, get_compression_type_indicated_by_filename( filename ))) ;
	}

	void RFormatStatSink::setup( std::auto_ptr< std::ostream > stream_ptr ) {
		set_stream( stream_ptr ) ;
	}
	
	void RFormatStatSink::write_value( double const& value ) {
		write_header_if_necessary() ;
		write_seperator_if_necessary() ;
		if( value == std::numeric_limits< double >::infinity() ) {
			stream() << "inf" ;
		} else {
			stream() << std::fixed << std::setprecision( m_precision ) << value ;
		}
	}
	
	void RFormatStatSink::write_descriptive_text() {
		stream() << m_descriptive_text ;
		if( m_descriptive_text.size() > 0 && m_descriptive_text[ m_descriptive_text.size() - 1] != '\n' ) {
			stream() << "\n" ;
		}
	}
	
	void RFormatStatSink::write_column_names() {
		for( std::size_t i = 0; i < number_of_columns(); ++i ) {
			stream() << column_names()[i] << " " ;
		}
	}
}
