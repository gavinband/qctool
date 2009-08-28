#ifndef NULL_OSTREAM_HPP
#define NULL_OSTREAM_HPP

#include <iostream>

struct null_ostream: public std::ostream
{
	struct null_streambuf: public std::streambuf
	{
		int overflow( int c ) { return std::char_traits< char >::not_eof(c) ; }
	} ;
	
	null_ostream(): std::ios( &m_buf), std::ostream( &m_buf ) {}

private:
	null_streambuf m_buf ;
} ;

#endif
