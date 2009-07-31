#include <vector>
#include <string>
#include <iostream>
#include <sstream>
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

SampleRow::SampleRow( SampleRow const& other ) :
	m_id1( other.m_id1 ),
	m_id2( other.m_id2 ),
	m_further_data( other.m_further_data ),
	m_column_headings( other.m_column_headings ),
	m_column_types( other.m_column_types )
{
}

void SampleRow::reset( std::vector<std::string> const& column_headings, std::vector<char> const& column_types ) {
	m_column_headings = column_headings ;
	m_column_types = column_types ;
	m_id1 = m_id2 = "" ;
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

void SampleRow::add_column( std::string const& heading, char type, double value ) {
	if( !have_column( heading )) {
		m_column_headings.push_back( heading ) ;
		m_column_types.push_back( type ) ;
		m_further_data[heading] = value ;
	}
}

bool SampleRow::have_column( std::string const& heading ) const {
	std::map< std::string, double >::const_iterator
		where = m_further_data.find( heading ) ;
	return ( where != m_further_data.end()) ;
}

void SampleRow::set_value( std::string const& heading, double value ) {
	assert( have_column( heading )) ;
	m_further_data[ heading ] = value ;
}

double SampleRow::get_double_value( std::string const& heading ) const {
	std::map< std::string, double > ::const_iterator
		where = m_further_data.find( heading ) ;
	assert( where != m_further_data.end() ) ;
	return where->second ;
}

std::string SampleRow::get_string_value( std::string const& name ) const {
	if( name == "ID1" ) {
		return m_id1 ;
	}
	else if( name == "ID2" ) {
		return m_id2 ;
	}
	else {
		double value = get_double_value( name ) ;
		std::ostringstream ostr ;
		ostr << value ;
		return ostr.str() ;
	}
}

std::istream& operator>>( std::istream& aStream, SampleRow& row ) {
	std::string line ;
	std::getline( aStream, line ) ;
	std::istringstream sStream( line ) ;
	sStream >> row.m_id1 ;
	sStream >> row.m_id2 ;

	std::size_t i = 2 ;
	while( i < row.m_column_headings.size() ) {
		double next_elt ;
		sStream >> next_elt ;
		row.m_further_data[ row.m_column_headings[i] ] = next_elt ;
		++i ;
	}

	if( sStream.fail() ) {
		aStream.setstate( std::ios::failbit ) ;
	}

	sStream.peek() ;
	if( !sStream.eof() ) {
		throw BadSampleRowFormatException( "Extra data at end of sample row" ) ;
	}

	return aStream ;
}

std::ostream& operator<<( std::ostream& aStream, SampleRow const& row ) {
	aStream << row.ID1() << " " << row.ID2() ;
	for( std::size_t i = 2 ; i < row.column_headings().size(); ++i ) {
		if( i > 0 )
			aStream << " " ;
		aStream << row.further_data( row.column_headings()[i] );
	}
	return aStream << "\n" ;
}
