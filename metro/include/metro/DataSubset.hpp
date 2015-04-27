
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_DATA_SUBSET_HPP
#define METRO_DATA_SUBSET_HPP

#include <iostream>
#include <vector>
#include "metro/DataRange.hpp"

namespace metro {
	
	// This class maintains a union of non-overlapping ranges of integers.
	// Used as a specification of ranges of data in likelihoods and elsewhere.
	struct DataSubset {
		DataSubset():
			m_size(0)
		{}

		DataSubset( DataRange const& range ):
			m_rep( 1, range ),
			m_size( range.size() )
		{}

		DataSubset( DataSubset const& other ):
			m_rep( other.m_rep ),
			m_size( other.m_size )
		{}

		DataSubset& operator=( DataSubset const& other ) {
			m_rep = other.m_rep ;
			m_size = other.m_size ;
			return *this ;
		}

		std::size_t size() const { return m_size ; }
		std::size_t number_of_subranges() const { return m_rep.size() ; }
		DataRange const& operator[]( std::size_t i ) const { return m_rep[i] ; }
		
		void add( DataRange const& range ) {
			add( DataSubset( range )) ;
		}

		void add( DataSubset const& other ) {
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
				i = j ;
			}
		}

	private:
		std::vector< DataRange > m_rep ;
		std::size_t m_size ;
	} ;
	
	bool operator== ( DataSubset const& left, DataSubset const& right ) ;
	bool operator< ( DataSubset const& left, DataSubset const& right ) ;
	
	std::ostream& operator<<( std::ostream& out, DataSubset const& range ) ;
}

#endif
