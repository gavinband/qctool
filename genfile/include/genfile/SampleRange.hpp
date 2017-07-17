
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SAMPLE_RANGE_HPP
#define GENFILE_SAMPLE_RANGE_HPP

#include <iostream>

namespace genfile {
	struct SampleRange {
		SampleRange() ;
		SampleRange( SampleRange const& other ) ;
		SampleRange( int begin, int end ) ;
		SampleRange& operator=( SampleRange const& other ) ;

		int begin() const { return m_begin ; }
		int end() const { return m_end ; }
		int size() const { return m_end - m_begin ; }

		private:
			int m_begin ;
			int m_end ;
	} ;
	
	bool operator== ( SampleRange const& left, SampleRange const& right ) ;
	bool operator< ( SampleRange const& left, SampleRange const& right ) ;
	
	std::ostream& operator<<( std::ostream& out, SampleRange const& range ) ;
}

#endif
