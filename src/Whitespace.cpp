#include <iostream>
#include "Whitespace.hpp"

Whitespace::Whitespace( std::string const whitespace_chars )
: m_whitespace_chars( whitespace_chars )
{}

std::istream& operator>>( std::istream& aStream, Whitespace & whitespace ) {
	char c;
	whitespace.m_buffer.clear() ;
	while( aStream.good() ) {
		c = aStream.peek() ;
		if( whitespace.m_whitespace_chars.find( c ) == std::string::npos ) {
			break ;
		}
		aStream.get() ;
		whitespace.m_buffer.append( 1, c ) ;
	}

	return aStream ;
}

