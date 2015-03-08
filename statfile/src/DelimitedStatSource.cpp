
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <limits>
#include "genfile/Error.hpp"
#include "genfile/string_utils/string_utils.hpp"
#include "statfile/DelimitedStatSource.hpp"

namespace statfile {
	DelimitedStatSource::DelimitedStatSource(
		std::auto_ptr< std::istream > stream_ptr,
		std::string delimiter,
		OptionalString const& ignore_until,
		OptionalString const& ignore_from
	): 	
		m_filename( "(unnamed stream)" ),
		m_comment_character( '#' ),
		m_delimiter( delimiter ),
		m_quotes( "\"\"" ),
		m_number_of_ignored_lines( 0 ),
		m_ignore_until( ignore_until ),
		m_ignore_from( ignore_from )
	{
		setup( stream_ptr ) ;
	}

	DelimitedStatSource::DelimitedStatSource(
		std::string const& filename,
		std::string delimiter,
		OptionalString const& ignore_until,
		OptionalString const& ignore_from
	): 	
		m_filename( filename ),
		m_comment_character( '#' ),
		m_delimiter( delimiter ),
		m_quotes( "\"\"" ),
		m_number_of_ignored_lines( 0 ),
		m_ignore_until( ignore_until ),
		m_ignore_from( ignore_from )
	{
		setup( filename ) ;
	}

	void DelimitedStatSource::reset_stream_to_start() {
		try {
			IstreamAggregator::reset_stream_to_start() ;
		}
		catch( std::ios_base::failure const& ) {
			set_stream( open_text_file_for_input( m_filename ) ) ;
			assert( stream() ) ;
		}
	}

	void DelimitedStatSource::reset_to_start() {
		reset_stream_to_start() ;
		// Skip comment and column header lines.
		std::string line ;
		for( std::size_t i = 0; i < m_number_of_comment_lines + m_number_of_ignored_lines + 1; ++i ) {
			std::getline( stream(), line ) ;
		}
		Base::reset_to_start() ;
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
		if( current_column() == 0 ) {
			read_one_line() ;
		}
	}

	void DelimitedStatSource::ignore_all() {
		if( current_column() == 0 ) {
			read_one_line() ;
		}
	}

	void DelimitedStatSource::read_one_line() {
		std::getline( stream(), m_current_line ) ;

		// handle comment lines by ignoring them.
		while( stream() && m_current_line.size() > 0 && m_current_line[0] == m_comment_character ) {
			std::getline( stream(), m_current_line ) ;
		}
			
		m_current_fields = split_line( m_current_line, m_delimiter, m_quotes ) ;
		if( m_current_fields.size() != number_of_columns() ) {
			throw genfile::MalformedInputError( get_source_spec(), number_of_rows_read() ) ;
		}
	}

	std::vector< genfile::string_utils::slice > DelimitedStatSource::split_line(
		std::string const& line,
		std::string const& delimiter,
		std::string const& quotes
	) const {
		std::vector< genfile::string_utils::slice > result ;
		int in_quote = 0 ;
		std::size_t last_i = 0 ;
		for( std::size_t i = 0; i < line.size(); ) {
			if( !in_quote && line.compare( i, delimiter.size(), delimiter ) == 0 ) {
				result.push_back( get_unquoted_substring( line, last_i, i - last_i, quotes )) ;
				i += delimiter.size() ;
				last_i = i ;
			}
			else if( line[i] == quotes[0] ) {
				in_quote = ( in_quote + 1 ) % 2 ;
				++i ;
			}
			else {
				++i ;
			}
		}
		if( in_quote ) {
			throw genfile::MalformedInputError( get_source_spec(), number_of_rows_read() ) ;
		}
		// handle trailing newline.
		if( line.size() > 0 && line[ line.size() - 1 ] == '\r' ) {
			result.push_back( get_unquoted_substring( line, last_i, line.size() - last_i - 1, quotes )) ;
		} else {
			result.push_back( get_unquoted_substring( line, last_i, line.size() - last_i, quotes )) ;
		}
		
		return result ;
	}
	
	genfile::string_utils::slice DelimitedStatSource::get_unquoted_substring(
		std::string const& big_string,
		std::size_t pos,
		std::size_t length,
		std::string const& quotes
	) {
		if( big_string[pos] == quotes[0] && big_string[pos] == quotes[1] ) {
			return genfile::string_utils::slice( big_string ).substr( pos + 1, pos + length - 1 ) ;
		}
		else {
			return genfile::string_utils::slice( big_string ).substr( pos, pos + length ) ;
		}
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
		m_number_of_ignored_lines = 0 ;
		if( m_ignore_until ) {
			std::string line ;
			bool ignoreUntilStringFound = false ;
			while( std::getline( stream(), line ) ) {
				++m_number_of_ignored_lines ;
				std::size_t end = line.size() ;
				if( line.size() > 0 && line[ line.size() -1 ] == '\r' ) {
					--end ;
				}
				if( line.compare( 0, end, m_ignore_until.get() ) == 0 ) {
					ignoreUntilStringFound = true ;
					break ;
				}
			}
			if( !ignoreUntilStringFound ) {
				throw genfile::MalformedInputError( 
					"DelimitedStatSource::setup()",
					"Expected line \"" + m_ignore_until.get() + "\" was not found",
					m_number_of_ignored_lines + m_number_of_comment_lines
				) ;
			}
		}
		turn_off_ios_exceptions() ;
		m_number_of_rows = count_remaining_lines() - 1 ;
		// Go back to the start and get the column names, prepare for reading.
		reset_stream_to_start() ;
		assert( stream() ) ;
		turn_on_ios_exceptions() ;
		m_descriptive_text = read_descriptive_comments() ;
		{
			std::string line ;
			for( std::size_t i = 0; i < m_number_of_ignored_lines; ++i ) {
				std::getline( stream(), line ) ;
			}
		}
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
				line = line.substr( 1, line.size() ) ; // skip comment char
			}
			if( line.size() > 0 && line[0] == ' ' ) {
				line = line.substr( 1, line.size() ) ; // skip space.
			}
			if( line.size() > 0 && line[ line.size() - 1 ] == '\r' ) {
				line = line.substr( 0, line.size() - 1 ) ; // remove trailing CR
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
		std::vector< genfile::string_utils::slice > column_names = split_line( line, m_delimiter, m_quotes ) ;
		for( std::size_t i = 0; i < column_names.size(); ++i ) {
			add_column( column_names[i] ) ;
		}
	}
	
	std::size_t DelimitedStatSource::count_remaining_lines() {
		std::string line ;
		std::size_t count = 0 ;
		bool foundIgnoreFromString = false ;
		
		while( std::getline( stream(), line )) {
			std::size_t end = line.size() ;
			if( line.size() > 0 && line[ end - 1 ] == '\r' ) {
				--end ;
			}
			if( m_ignore_from && line.compare( 0, end, m_ignore_from.get() ) == 0 ) {
				foundIgnoreFromString = true ;
				break ;
			} else if( line.size() == 0 || line[0] != m_comment_character ) {
				++count ;
			}
		}
		if( m_ignore_from && !foundIgnoreFromString ) {
			throw genfile::MalformedInputError( 
				"DelimitedStatSource::setup()",
				"Expected line \"" + m_ignore_from.get() + "\" was not found",
				m_number_of_ignored_lines + m_number_of_comment_lines + count
			) ;
		}
		return count ;
	}
	
	std::string DelimitedStatSource::get_source_spec() const {
		return m_filename ;
	}
	
	// Specialisation of do_read_value for strings, since we don't need the istringstream here.
	template<>
	void DelimitedStatSource::do_read_value< std::string >( std::string& value ) {
		//std::cerr << "do_read_value: size of line is " << m_current_fields.size() << ".\n" ;
		if( current_column() == 0 ) {
			read_one_line() ;
		}
		std::string baked = m_current_fields[ current_column() ] ;
	 	value.swap( baked );
	}

	// Specialisation of do_read_value for doubles, to deal with infinities
	// represented in the file as "inf".
	template<>
	void DelimitedStatSource::do_read_value< double >( double& value ) {
		if( current_column() == 0 ) {
			read_one_line() ;
		}
		genfile::string_utils::slice const& elt = m_current_fields[ current_column() ] ;
		if( elt.size() == 2 && elt[0] == 'N' && elt[1] == 'A' ) {
			value = std::numeric_limits< double >::quiet_NaN() ;
		}
		else {
			value = genfile::string_utils::to_repr< double >( elt ) ;
		}
	}
}
