
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef __GTOOL_WHITESPACE__
#define __GTOOL_WHITESPACE__

#include <iostream>
#include <string>

struct Whitespace
{
	Whitespace( std::string const whitespace_chars = " \t\n\r" ) ;

	std::size_t size() const { return m_buffer.size() ; }
	bool empty() const { return m_buffer.empty() ; }

	private:

		std::string m_whitespace_chars ;
		std::string m_buffer ;

	// Read whitespace from the input stream.
	// Only consumes what is read.
	friend std::istream& operator>>( std::istream&, Whitespace & ) ;
} ;

#endif

