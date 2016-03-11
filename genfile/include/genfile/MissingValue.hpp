
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_MISSINGVALUE_HPP
#define GENFILE_MISSINGVALUE_HPP

#include <iostream>

namespace genfile {
	struct MissingValue
	{
		MissingValue() {}
		MissingValue( MissingValue const& other ) {}
		MissingValue& operator=( MissingValue const& other ) { return *this ; }
		bool operator<( MissingValue const& other ) const ;
		bool operator<=( MissingValue const& other ) const ;
		bool operator==( MissingValue const& other ) const ;
	} ;

	std::ostream& operator<<( std::ostream& o, MissingValue const& v ) ;
}

#endif
