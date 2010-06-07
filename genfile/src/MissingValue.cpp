#include <iostream>
#include "genfile/MissingValue.hpp"

namespace genfile {
	std::ostream& operator<<( std::ostream& o, MissingValue const& v ) {
		return o << "NA" ;
	}
	
	bool MissingValue::operator==( MissingValue const& other ) const {
		return true ;
	}

	bool MissingValue::operator<( MissingValue const& other ) const {
		return false ;
	}
}
