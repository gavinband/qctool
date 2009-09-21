#ifndef SET_OF_POSITIONS_HPP
#define SET_OF_POSITIONS_HPP

#include "Position.hpp"

struct SetOfPositions
{
	virtual ~SetOfPositions() {} ;
	virtual bool contains( Position const position ) const = 0 ;
} ;

#endif
