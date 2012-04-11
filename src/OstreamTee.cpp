
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <map>
#include <string>
#include <cassert>
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
