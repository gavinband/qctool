#ifndef GENFILE_MISSINGVALUE_HPP
#define GENFILE_MISSINGVALUE_HPP

#include <iostream>

namespace genfile {
	struct MissingValue
	{
		bool operator<( MissingValue const& other ) const ;
		bool operator==( MissingValue const& other ) const ;
	} ;

	std::ostream& operator<<( std::ostream& o, MissingValue const& v ) ;
}

#endif
