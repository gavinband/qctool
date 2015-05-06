//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <vector>
#include "metro/DataRange.hpp"
#include "metro/DataSubset.hpp"

namespace metro {

	void DataSubset::add( DataSubset const& other ) {
		m_rep.insert( m_rep.end(), other.m_rep.begin(), other.m_rep.end() ) ;
		std::sort( m_rep.begin(), m_rep.end() ) ;
		// Now join ranges that need to be joined.
		m_size = 0 ;
		for( std::size_t i = 0; i < m_rep.size(); ++i ) {
			int begin = m_rep[i].begin() ;
			int end = m_rep[i].end() ;
			std::size_t j = i ;
			for( ++j; j < m_rep.size() && m_rep[j].begin() <= end; ++j ) {
				end = std::max( end, m_rep[j].end() ) ;
			}
			m_rep[i] = DataRange( begin, end ) ;
			m_rep.erase( m_rep.begin() + i + 1, m_rep.begin() + j ) ;
			m_size += m_rep[i].size() ;
		}
	}

	std::ostream& operator<<( std::ostream& out, DataSubset const& range ) {
		for( std::size_t i = 0; i < range.number_of_subranges(); ++i ) {
			out << ( i > 0 ? "," : "" ) << range[i] ;
		}
		return out ;
	}
}
