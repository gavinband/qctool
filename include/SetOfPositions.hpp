
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SET_OF_POSITIONS_HPP
#define SET_OF_POSITIONS_HPP

#include "Position.hpp"

struct SetOfPositions
{
	virtual ~SetOfPositions() {} ;
	virtual bool contains( Position const position ) const = 0 ;
} ;

#endif
