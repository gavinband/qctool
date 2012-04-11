
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CASE_CONTROL_TOOLS_POSITION_RANGE_HPP
#define CASE_CONTROL_TOOLS_POSITION_RANGE_HPP

#include <iostream>
#include <stddef.h>
#include "Position.hpp"
#include "SetOfPositions.hpp"

struct PositionRange: public SetOfPositions
{
	// empty range.
	PositionRange()
		: m_empty( true )
	{}
	
	PositionRange( Position start, Position end )
		: m_start( start ),
		  m_end( end ),
		  m_empty( false )
	{
		assert( m_start <= m_end ) ;
	}

	PositionRange( PositionRange const& other )
		: m_start( other.m_start ),
		  m_end( other.m_end ),
		  m_empty( other.m_empty )
	{}

	Position const& start() const { assert( !is_empty() ) ; return m_start ;}
	Position const& end() const { assert( !is_empty() ) ; return m_end ;}
	std::size_t size() const {
		if( is_empty()) {
			return 0u ;
		}
		else {
			return m_end + 1 - m_start ;
		}
	}

	PositionRange with_margin( Position margin ) const {
		assert( !is_empty() ) ;
		return PositionRange( ((m_start < margin) ? 0 : m_start - margin), m_end + margin ) ;
	}
	
	bool is_empty() const { return m_empty ; }

	bool contains( Position const position ) const {
		return (!m_empty) && (position >= m_start) && (position <= m_end) ;
	}

private:
	Position m_start ;
	Position m_end ;
	bool m_empty ;	
} ;

bool operator<( PositionRange const& left, PositionRange const& right ) {
	return left.start() < right.start() ;
}

std::ostream& operator<<( std::ostream& out, PositionRange const& range ) {
	return out << "[" << range.start() << "," << range.end() << "]" ;
}

typedef std::vector< PositionRange > PositionRangeList ;
std::ostream& operator<<( std::ostream& out, PositionRangeList const& ranges ) {
	for( std::size_t i = 0; i < ranges.size(); ++i ) {
		if( i > 0 ) {
			out << "," ;
		}
		out << ranges[i] ;
	}
	return out ;
}

PositionRange intersection( PositionRange const& left, PositionRange const& right ) {
	Position left_end = std::max( left.start(), right.start() ) ;
	Position right_end = std::min( left.end(), right.end() ) ;
	if( left_end > right_end ) {
		return PositionRange() ;
	}
	else {
		return PositionRange( left_end, right_end ) ;
	}
}


#endif
