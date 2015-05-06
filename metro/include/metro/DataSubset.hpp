
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

		void add( DataSubset const& other ) ;

	private:
		std::vector< DataRange > m_rep ;
		std::size_t m_size ;
	} ;
	
	std::ostream& operator<<( std::ostream& out, DataSubset const& range ) ;
}

#endif
