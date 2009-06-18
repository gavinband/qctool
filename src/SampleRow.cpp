#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include "SampleRow.hpp"
#include "string_utils.hpp"
#include "Whitespace.hpp"

bool check_char_is_space( char c ) { return c == ' ' ; }
bool check_char_is_newline( char c ) { return c == '\n' ; }

SampleRow::SampleRow() {}

SampleRow::SampleRow( std::vector< std::string > const& column_headings, std::vector<char> const& column_types )
{
	reset( column_headings, column_types ) ;
}

void SampleRow::reset( std::vector<std::string> const& column_headings, std::vector<char> const& column_types ) {
	m_column_headings = column_headings ;
	m_column_types = column_types ;
	m_id1 = m_id2 = "" ;
	m_missing = 0.0 ;
	m_further_data.clear() ;

	assert( m_column_headings.size() >= 3 ) ;
	assert( m_column_types.size() == m_column_headings.size()) ;
	assert( m_column_headings[2] == "missing" ) ;
	assert( m_column_types[0] == '0' ) ;
	assert( m_column_types[1] == '0' ) ;
	assert( m_column_types[2] == '0' ) ;
}

double SampleRow::further_data( std::string const& column_heading ) const {
	std::map< std::string, double >::const_iterator
		further_data_i = m_further_data.find( column_heading ) ;
	assert( further_data_i != m_further_data.end() ) ;
	return further_data_i->second ;
}

std::istream& operator>>( std::istream& aStream, SampleRow& row ) {
	aStream >> std::noskipws ;
	char c ;
	
	aStream >> row.m_id1 >> c ;
	if( !check_char_is_space(c) ) {
		throw BadSampleRowFormatException( "Expected single space after sample row id 1") ;
	};
	aStream >> row.m_id2 >> c ;
	if( !check_char_is_space(c) ) {
		throw BadSampleRowFormatException( "Expected single space after sample row id 2") ;
	};
	aStream >> row.m_missing ;

	std::size_t i = 3 ;
	while( i < row.m_column_headings.size() ) {
		aStream >> c ;
		if( !check_char_is_space(c) ) {
			throw BadSampleRowFormatException( "Expected single space before sample row entry.") ;
		};
		double next_elt ;
		aStream >> next_elt ;
		row.m_further_data[ row.m_column_headings[i] ] = next_elt ;
		++i ;
	}

	// At the end of the row, allow trailing whitespace.
	// This should be followed by a newline or the eof.
	Whitespace whitespace( " \t" ) ;
	aStream >> whitespace ;
	aStream.peek() ; // ensure eof is flagged if we've reached the end of the file.
	if( !aStream.eof() ) {
		if( check_char_is_newline( static_cast<char>( aStream.peek() ))) {
			aStream.get() ;
			aStream.peek() ; // ensure eof is flagged if we've reached the end of the file.
		}
		else {
			std::string msg = "Expected newline or eof at end of sample row, got \'" ;
			msg += static_cast< char >( aStream.peek() ) ;
			msg += "\'." ;
			throw BadSampleRowFormatException( msg ) ;
		}
	}
	
	return aStream >> std::skipws ;
}

std::ostream& operator<<( std::ostream& aStream, SampleRow const& row ) {
	aStream << row.ID1() << " " << row.ID2() << " " << row.missing() ;
	for( std::size_t i = 3 ; i < row.column_headings().size(); ++i ) {
		if( i > 0 )
			aStream << " " ;
		aStream << row.further_data( row.column_headings()[i] );
	}
	return aStream << "\n" ;
}
