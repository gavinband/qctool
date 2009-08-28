#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "statfile/RFormatStatSource.hpp"

namespace statfile {
	RFormatStatSource::RFormatStatSource( std::auto_ptr< std::istream > stream_ptr )
	: m_comment_character( '#' )
	{
		setup( stream_ptr ) ;
	}

	RFormatStatSource::RFormatStatSource( std::string const& filename )
	: m_comment_character( '#' )
	{
		setup( filename ) ;
	}

	void RFormatStatSource::read_value( int32_t& field ) {
		stream() >> field ;
	}

	void RFormatStatSource::read_value( uint32_t& field ) {
		stream() >> field ;
	}

	void RFormatStatSource::read_value( std::string& field ) {
		stream() >> field ;
	}

	void RFormatStatSource::read_value( double& field ) {
		std::string str_field ;
		stream() >> str_field ;
		if( str_field == "inf" ) {
			field = std::numeric_limits< double >::infinity() ;
		} else {
			std::istringstream istr( str_field ) ;
			istr >> field ;
			istr.peek() ;
			assert( istr.eof()) ;	
		}
	}

	void RFormatStatSource::ignore_value() {
		std::string field ;
		stream() >> field ;
	}

	void RFormatStatSource::ignore_all() {
		std::string line ;
		std::getline( stream(), line ) ;
	}

	void RFormatStatSource::setup( std::auto_ptr< std::istream > stream_ptr ) {
		// First count the lines in the file.
		set_stream( stream_ptr ) ;
		read_descriptive_comments() ;
		m_total_number_of_rows = count_remaining_lines() - 1 ;
		// Go back to the start and get the column names, prepare for reading.
		reset_stream_to_start() ;
		assert( stream() ) ;
		m_descriptive_text = read_descriptive_comments() ;
		read_column_names() ;
	}

	void RFormatStatSource::setup( std::string const& filename ) {
		setup( open_text_file_for_input( filename, get_compression_type_indicated_by_filename( filename ))) ;
	}

	std::string RFormatStatSource::read_descriptive_comments() {
		std::string line, result ;
		while( stream().peek() == m_comment_character || stream().peek() == '\n' ) {
			std::getline( stream(), line ) ;
			if( line.size() > 0 ) {
				line = line.substr( 1, line.size() ) ; // miss off comment char.
			}
			if( result.size() > 0 ) {
				result += '\n' ;
			}
			result += line ;
		}
		return result ;
	}

	void RFormatStatSource::read_column_names() {
		std::string line ;
		std::getline( stream(), line ) ;
		std::istringstream istr( line ) ;
		std::string elt ;
		while( istr >> elt ) {
			add_column( elt ) ;
		}
	}
	
	std::size_t RFormatStatSource::count_remaining_lines() {
		std::string line ;
		std::size_t count = 0 ;
		while( std::getline( stream(), line )) {
			++count ;
		}
		return count ;
	}
}
