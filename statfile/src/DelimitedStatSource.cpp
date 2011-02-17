#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <limits>

#include "statfile/DelimitedStatSource.hpp"

namespace statfile {
	DelimitedStatSource::DelimitedStatSource(
		std::auto_ptr< std::istream > stream_ptr,
		std::string delimiter
	)
	: 	m_comment_character( '#' ),
		m_delimiter( delimiter ),
		m_strip_chars( "\"" )
	{
		setup( stream_ptr ) ;
	}

	DelimitedStatSource::DelimitedStatSource( std::string const& filename, std::string delimiter )
	: 	m_comment_character( '#' ),
		m_delimiter( delimiter ),
		m_strip_chars( "\"" )
	{
		setup( filename ) ;
	}

	void DelimitedStatSource::reset_to_start() {
		reset_stream_to_start() ;
		// Skip comment and column header lines.
		std::string line ;
		for( std::size_t i = 0; i < m_number_of_comment_lines + 1; ++i ) {
			std::getline( stream(), line ) ;
		}
		base_t::reset_to_start() ;
	}

	void DelimitedStatSource::read_value( int32_t& field ) {
		do_read_value( field ) ;
	}

	void DelimitedStatSource::read_value( uint32_t& field ) {
		do_read_value( field ) ;
	}

	void DelimitedStatSource::read_value( std::string& field ) {
		do_read_value( field ) ;
	}

	void DelimitedStatSource::read_value( double& field ) {
		do_read_value( field ) ;	
	}

	void DelimitedStatSource::ignore_value() {
		// nothing to do.
	}

	void DelimitedStatSource::ignore_all() {
		// nothing to do.
	}

	void DelimitedStatSource::read_one_line() {
		std::string line ;
		std::getline( stream(), line ) ;
		m_current_fields = split_line( line, m_delimiter, m_strip_chars ) ;
		if( m_current_fields.size() != number_of_columns() ) {
			throw FileIsInvalidError() ;
		}
	}

	std::vector< std::string > DelimitedStatSource::split_line(
		std::string const& line,
		std::string const& delimiter,
		std::string const& strip_chars
	) {
		std::vector< std::string > result ;
		std::size_t begin_pos = 0, delim_pos ;
		do {
			delim_pos = line.find( delimiter, begin_pos ) ;
			if( delim_pos == std::string::npos )
				delim_pos = line.size() ;
			result.push_back( genfile::string_utils::strip( line.substr( begin_pos, delim_pos - begin_pos ), strip_chars )) ;
			begin_pos = delim_pos + delimiter.size() ;
		}
		while( delim_pos != line.size() ) ;
		return result ;
	}
	
	std::string DelimitedStatSource::strip( std::string const& string_to_strip, std::string const& strip_chars ) {
		std::size_t lpos = string_to_strip.find_first_not_of( strip_chars ) ;
		return ( lpos == std::string::npos )
			? ""
			: string_to_strip.substr( lpos, string_to_strip.find_last_not_of( strip_chars ) - lpos + 1 ) ;
	}

	void DelimitedStatSource::setup( std::auto_ptr< std::istream > stream_ptr ) {
		// First count the lines in the file.
		set_stream( stream_ptr ) ;
		read_descriptive_comments() ;
		turn_off_ios_exceptions() ;
		m_number_of_rows = count_remaining_lines() - 1 ;
		// Go back to the start and get the column names, prepare for reading.
		reset_stream_to_start() ;
		assert( stream() ) ;
		turn_on_ios_exceptions() ;
		m_descriptive_text = read_descriptive_comments() ;
		read_column_names() ;
	}

	void DelimitedStatSource::setup( std::string const& filename ) {
		setup( open_text_file_for_input( filename, get_compression_type_indicated_by_filename( filename ))) ;
	}

	std::string DelimitedStatSource::read_descriptive_comments() {
		std::string line, result ;
		m_number_of_comment_lines = 0 ;
		while( stream().peek() == m_comment_character || stream().peek() == '\n' ) {
			std::getline( stream(), line ) ;
			if( line.size() > 0 ) {
				line = line.substr( 1, line.size() ) ; // miss off comment char.
			}
			if( result.size() > 0 ) {
				result += '\n' ;
			}
			result += line ;
			++m_number_of_comment_lines ;
		}
		return result ;
	}

	void DelimitedStatSource::read_column_names() {
		std::string line ;
		std::getline( stream(), line ) ;
		std::vector< std::string > column_names = split_line( line, m_delimiter, m_strip_chars ) ;
		for( std::size_t i = 0; i < column_names.size(); ++i ) {
			add_column( column_names[i] ) ;
		}
	}
	
	std::size_t DelimitedStatSource::count_remaining_lines() {
		std::string line ;
		std::size_t count = 0 ;
		while( std::getline( stream(), line )) {
			++count ;
		}
		return count ;
	}
	
	// Specialisation of do_read_value for strings, since we don't need the istringstream here.
	template<>
	void DelimitedStatSource::do_read_value< std::string >( std::string& value ) {
		if( current_column() == 0 ) {
			read_one_line() ;
		}
	 	value = m_current_fields[ current_column() ] ;
	}

	// Specialisation of do_read_value for doubles, to deal with infinities
	// represented in the file as "inf".
	template<>
	void DelimitedStatSource::do_read_value< double >( double& value ) {
		std::string str_field ;
		do_read_value< std::string >( str_field ) ;
		if( str_field == "inf" ) {
			value = std::numeric_limits< double >::infinity() ;
		} else {
			std::istringstream istr( str_field ) ;
			istr >> value ;
			istr.peek() ;
			assert( istr.eof()) ;	
		}
	}
	
}
