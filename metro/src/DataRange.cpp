
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include <iostream>
#include "metro/DataRange.hpp"

namespace metro {
	DataRange::DataRange():
		m_begin( 0 ),
		m_end( 0 )
	{}

	DataRange::DataRange( DataRange const& other ):
		m_begin( other.m_begin ),
		m_end( other.m_end )
	{}

	DataRange::DataRange( int begin, int end ):
		m_begin( begin ),
		m_end( end )
	{
		assert( m_end >= m_begin ) ;
	}

	DataRange& DataRange::operator=( DataRange const& other ) {
		m_begin = other.m_begin ;
		m_end = other.m_end ;
		return *this ;
	}

	bool operator== ( DataRange const& left, DataRange const& right ) {
		return right.begin() == left.begin() && right.end() == left.end() ;
	}
	
	bool operator< ( DataRange const& left, DataRange const& right ) {
		return
			( left.begin() < right.begin() )
			| 
			( left.begin() == right.begin() && left.end() < right.end() )
		;
	}
	
	std::ostream& operator<<( std::ostream& out, DataRange const& range ) {
		return out << "[" << range.begin() << "-" << range.end() << ")" ;
	}
}
