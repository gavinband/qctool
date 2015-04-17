
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_DATA_RANGE_HPP
#define METRO_DATA_RANGE_HPP

#include <iostream>

namespace metro {
	struct DataRange {
		DataRange() ;
		DataRange( DataRange const& other ) ;
		DataRange( int begin, int end ) ;
		DataRange& operator=( DataRange const& other ) ;

		int begin() const { return m_begin ; }
		int end() const { return m_end ; }
		int size() const { return m_end - m_begin ; }

		private:
			int m_begin ;
			int m_end ;
	} ;
	
	bool operator== ( DataRange const& left, DataRange const& right ) ;
	bool operator< ( DataRange const& left, DataRange const& right ) ;
	
	std::ostream& operator<<( std::ostream& out, DataRange const& range ) ;
}

#endif
