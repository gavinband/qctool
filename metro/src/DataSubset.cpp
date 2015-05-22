//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <vector>
#include <algorithm>
#include "metro/DataRange.hpp"
#include "metro/DataSubset.hpp"

namespace metro {

#if 0
	void DataSubset::add( DataRange const& other ) {
		std::vector< DataRange > newRep ;
		Rep::const_iterator i = m_rep.begin() ;
		for(
			;
			i != m_rep.end() && i->end() < other.begin();
			++i
		) ;
		std::copy( m_rep.begin(), i, std::back_inserter< newRep >() ) ;
		int begin = other.begin() ;
		int end = other.end() ;
		for(
			;
			i != m_rep.end() && i->begin() < other.end();
			++i
		) {
			begin = std::min( begin, i->begin() ) ;
			end = std::max( end, i->end() ) ;
		}



		for( std::size_t i = 0; i < m_rep.size() && ; ++i ) {

			if( m_rep[i].end() < m_rep[i].begin() || m_rep[i].begin() > m_rep[i].end() ) {
				end = i ;
			} else {
				// join
			}
		}
	}
#endif
	
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
