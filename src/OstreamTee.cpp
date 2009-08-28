#include <iostream>
#include <map>
#include <string>
#include "OstreamTee.hpp"

void OstreamTee::add_stream( std::string const& name, std::ostream& stream ) {
	m_streams[ name ] = &stream ;
}

std::ostream& OstreamTee::operator[]( std::string const& name ) {
	std::map< std::string, std::ostream* >::iterator
		where = m_streams.find( name ) ;
	assert( where != m_streams.end() ) ;
	return *(where->second) ;
}
