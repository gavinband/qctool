
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include <iostream>
#include "metro/SampleRange.hpp"

namespace metro {
	SampleRange::SampleRange():
		m_begin( 0 ),
		m_end( 0 )
	{}

	SampleRange::SampleRange( SampleRange const& other ):
		m_begin( other.m_begin ),
		m_end( other.m_end )
	{}

	SampleRange::SampleRange( int begin, int end ):
		m_begin( begin ),
		m_end( end )
	{
		assert( m_end >= m_begin ) ;
	}

	SampleRange& SampleRange::operator=( SampleRange const& other ) {
		m_begin = other.m_begin ;
		m_end = other.m_end ;
		return *this ;
	}

	bool operator== ( SampleRange const& left, SampleRange const& right ) {
		return right.begin() == left.begin() && right.end() == left.end() ;
	}
	
	bool operator< ( SampleRange const& left, SampleRange const& right ) {
		return
			( left.begin() < right.begin() )
			| 
			( left.begin() == right.begin() && left.end() < right.end() )
		;
	}
	
	std::ostream& operator<<( std::ostream& out, SampleRange const& range ) {
		return out << "[" << range.begin() << "-" << range.end() << ")" ;
	}

	std::ostream& operator<<( std::ostream& out, std::vector< SampleRange > const& ranges ) {
		for( std::size_t i = 0; i < ranges.size(); ++i ) {
			out << (i>0?",":"" ) << ranges[i] ;
		}
		return out ;
	}
}
