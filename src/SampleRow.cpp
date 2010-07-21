#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include "SampleRow.hpp"
#include "string_utils/string_utils.hpp"
#include "Whitespace.hpp"
#include "genfile/CohortIndividualSource.hpp"

bool check_char_is_space( char c ) { return c == ' ' ; }
bool check_char_is_newline( char c ) { return c == '\n' ; }

SampleRow::SampleRow() {
	m_column_headings.push_back( "id_1" ) ;
	m_column_headings.push_back( "id_2" ) ;
	m_column_headings.push_back( "missing" ) ;
	m_column_types.push_back( '0' ) ;
	m_column_types.push_back( '0' ) ;
	m_column_types.push_back( '0' ) ;
	m_further_data[ "id_1" ] = "NA" ;
	m_further_data[ "id_2" ] = "NA" ;
	m_further_data[ "missing" ] = "NA" ;
}

SampleRow::SampleRow( SampleRow const& other ) :
	m_further_data( other.m_further_data ),
	m_column_headings( other.m_column_headings ),
	m_column_types( other.m_column_types )
{
}

void SampleRow::reset( std::vector< genfile::CohortIndividualSource::SingleColumnSpec > const& column_spec ) {
	m_column_headings.resize( column_spec.size() ) ;
	m_column_types.resize( column_spec.size() ) ;
	for( std::size_t i = 0; i < column_spec.size(); ++i ) {
		m_column_headings[i] = column_spec[i].name() ;
		std::ostringstream ostr ;
		ostr << column_spec[i].type() ;
		assert( ostr.str().size() == 1 ) ;
		m_column_types[i] = ostr.str()[0] ;
	}
	m_further_data.clear() ;

	assert( m_column_headings.size() >= 3 ) ;
	assert( m_column_types.size() == m_column_headings.size()) ;
	assert( m_column_headings[2] == "missing" ) ;
	assert( m_column_types[0] == '0' ) ;
	assert( m_column_types[1] == '0' ) ;
	assert( m_column_types[2] == '0' ) ;
}

std::string SampleRow::ID1() const {
	FurtherData::const_iterator where = m_further_data.find( m_column_headings[0] ) ;
	assert( where != m_further_data.end() ) ;
	return where->second.as< std::string >() ;
}
std::string SampleRow::ID2() const {
	FurtherData::const_iterator where = m_further_data.find( m_column_headings[1] ) ;
	assert( where != m_further_data.end() ) ;
	return where->second.as< std::string >() ;
}

SampleRow::Entry SampleRow::further_data( std::string const& column_heading ) const {
	FurtherData::const_iterator further_data_i = m_further_data.find( column_heading ) ;
	assert( further_data_i != m_further_data.end() ) ;
	return further_data_i->second ;
}

void SampleRow::add_column( std::string const& heading, char type, Entry value ) {
	if( !has_column( heading )) {
		// A kludge.
		// To avoid breaking the order of rows (covariates first)
		// We insert all new columns as the 4th column.
		assert( m_column_headings.size() >= 3 ) ;
		assert( m_column_types.size() >= 3 ) ;
		std::vector<std::string>::iterator fourth_column_heading = m_column_headings.begin() + 3 ;
		std::vector<char>::iterator fourth_column_type = m_column_types.begin() + 3 ;
		m_column_headings.insert( fourth_column_heading, heading ) ;
		m_column_types.insert( fourth_column_type, type ) ;
		m_further_data[heading] = value ;
	}
}

bool SampleRow::has_column( std::string const& heading ) const {
	std::vector< std::string >::const_iterator where = std::find( m_column_headings.begin(), m_column_headings.end(), heading ) ;
	return ( where != m_column_headings.end()) ;
}

bool SampleRow::has_value( std::string const& name ) const {
	return SampleRow::has_column( name ) ;
}

void SampleRow::set_value( std::string const& heading, Entry value ) {
	assert( has_column( heading )) ;
	m_further_data[ heading ] = value ;
}

double SampleRow::get_double_value( std::string const& heading ) const {
	FurtherData::const_iterator where = m_further_data.find( heading ) ;
	assert( where != m_further_data.end() ) ;
	return where->second.as< double >() ;
}

std::string SampleRow::get_string_value( std::string const& heading ) const {
	std::ostringstream ostr ;
	FurtherData::const_iterator where = m_further_data.find( heading ) ;
	assert( where != m_further_data.end() ) ;
	ostr << where->second ;
	return ostr.str() ;
}

void SampleRow::read_ith_sample_from_source( std::size_t sample_i, genfile::CohortIndividualSource const& source ) {
	assert( sample_i < source.get_number_of_individuals() ) ;
	std::vector< genfile::CohortIndividualSource::SingleColumnSpec > column_spec = source.get_column_spec() ;
	reset( column_spec ) ;
	for( std::size_t i = 0; i < source.get_number_of_columns(); ++i ) {
		set_value( column_spec[i].name(), source.get_entry( sample_i, column_spec[i].name() )) ;
	}
} 

std::ostream& operator<<( std::ostream& aStream, SampleRow const& row ) {
	for( std::size_t i = 0 ; i < row.column_headings().size(); ++i ) {
		if( i > 0 )
			aStream << " " ;
		aStream << row.further_data( row.column_headings()[i] );
	}
	return aStream << "\n" ;
}
