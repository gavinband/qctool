
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GEN_TOOLS_OSTREAM_TEE_HPP
#define GEN_TOOLS_OSTREAM_TEE_HPP

#include <iostream>
#include <map>
#include <string>

class OstreamTee: public std::ostream  {
public:

	void add_stream( std::string const& name, std::ostream& stream ) ; 
	std::ostream& operator[]( std::string const& name ) ;

	template< typename T >
	friend OstreamTee& operator<<( OstreamTee & ostream_tee, T const& t ) ;

private:
	std::map< std::string, std::ostream* > m_streams ;
} ;

template< typename T >
OstreamTee& operator<<( OstreamTee& ostream_tee, T const& t ) {
	for(
		std::map< std::string, std::ostream* >::const_iterator i = ostream_tee.m_streams.begin() ;
		i != ostream_tee.m_streams.end() ;
		++i
	) {
		(*(i->second)) << t ;
	}

	return ostream_tee ;
}

#endif