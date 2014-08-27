
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <memory>
#include <limits>
#include "statfile/DelimitedStatSink.hpp"

namespace statfile {
	DelimitedStatSink::DelimitedStatSink( std::auto_ptr< std::ostream > stream_ptr, std::string const& delimiter )
	: 	m_comment_character( '#' ),
		m_delimiter( delimiter ),
		have_written_header( false ),
		m_precision(15)
	{
		setup( stream_ptr ) ;
	}

	DelimitedStatSink::DelimitedStatSink( std::string const& filename, std::string const& delimiter )
	: 	m_comment_character( '#' ),
		m_delimiter( delimiter ),
		have_written_header( false ),
		m_precision(15)
	{
		setup( filename ) ;
	}

	void DelimitedStatSink::set_descriptive_text( std::string const& text ) {
		m_descriptive_text = "# " ;
		for( std::string::const_iterator i = text.begin(); i != text.end(); ++i ) {
			m_descriptive_text += *i ;
			if( *i == '\n' ) {
				m_descriptive_text.append( "# " ) ;
			}
		}
	}

	void DelimitedStatSink::setup( std::string const& filename ) {
		setup( open_text_file_for_output( filename, get_compression_type_indicated_by_filename( filename ))) ;
	}

	void DelimitedStatSink::setup( std::auto_ptr< std::ostream > stream_ptr ) {
		set_stream( stream_ptr ) ;
	}
	
	void DelimitedStatSink::write_value( double const& value ) {
		write_header_if_necessary() ;
		write_seperator_if_necessary() ;
		if( value == std::numeric_limits< double >::infinity() ) {
			stream() << "inf" ;
		} else if( value == value ) {
			stream() << std::setprecision( m_precision ) << value ;
		} else {
			stream() << "NA" ;
		}
	}

	void DelimitedStatSink::write_descriptive_text() {
		stream() << m_descriptive_text ;
		if( m_descriptive_text.size() > 0 && m_descriptive_text[ m_descriptive_text.size() - 1] != '\n' ) {
			stream() << "\n" ;
		}
	}
	
	void DelimitedStatSink::write_column_names() {
		for( std::size_t i = 0; i < number_of_columns(); ++i ) {
			if( i > 0 ) {
				stream() << m_delimiter ;
			}
			stream() << column_names()[i] ;
		}
	}
}
