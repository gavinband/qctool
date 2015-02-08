
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
		m_always_escape_strings( false ),
		m_precision(6)
	{
		assert( m_delimiter.size() == 1 ) ;
		setup( stream_ptr ) ;
	}

	DelimitedStatSink::DelimitedStatSink( std::string const& filename, std::string const& delimiter )
	: 	m_comment_character( '#' ),
		m_delimiter( delimiter ),
		m_always_escape_strings( false ),
		m_precision(6)
	{
		assert( m_delimiter.size() == 1 ) ;
		setup( filename ) ;
	}

	void DelimitedStatSink::set_descriptive_text( std::string const& text ) {
		m_descriptive_text = std::string( 1, m_comment_character ) + " " ;
		for( std::string::const_iterator i = text.begin(); i != text.end(); ++i ) {
			m_descriptive_text += *i ;
			if( *i == '\n' ) {
				m_descriptive_text.append( std::string( 1, m_comment_character ) + " " ) ;
			}
		}
	}

	void DelimitedStatSink::setup( std::string const& filename ) {
		setup( open_text_file_for_output( filename, get_compression_type_indicated_by_filename( filename ))) ;
	}

	void DelimitedStatSink::setup( std::auto_ptr< std::ostream > stream_ptr ) {
		set_stream( stream_ptr ) ;
	}

	void DelimitedStatSink::begin_data_impl() {
		write_header_if_necessary() ;
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
	
	namespace {
		std::string escape_string(
			std::string const& value,
			char const delimiter,
			std::string const& begin_quote,
			std::string const& end_quote,
			bool always_escape = false
		) {
			std::string result ;
			if( always_escape || ( value.find_first_of( delimiter ) != std::string::npos ) ) {
				assert( value.find_first_of( begin_quote ) == std::string::npos ) ; // embedded quotes not supported right now.
				result = begin_quote + value + end_quote ;
			} else {
				result = value ;
			}
			return result ;
		}
	}
	void DelimitedStatSink::write_value( std::string const& value ) {
		// if delimiter is a single char, we escape it.
		write_header_if_necessary() ;
		write_seperator_if_necessary() ;
		std::string output_value = escape_string( value, m_delimiter[0], "\"", "\"", m_always_escape_strings ) ;
		stream() << output_value ;
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
			stream() << escape_string( column_names()[i], m_delimiter[0], "\"", "\"", m_always_escape_strings ) ;
		}
	}
}
