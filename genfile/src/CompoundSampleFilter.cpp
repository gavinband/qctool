
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include <memory>
#include <vector>
#include <set>
#include <boost/ptr_container/ptr_vector.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/SampleFilter.hpp"
#include "genfile/CompoundSampleFilter.hpp"

namespace genfile {
	void CompoundSampleFilter::add_clause( SampleFilter::UniquePtr clause ) {
		m_clauses.push_back( clause ) ;
	}

	std::size_t CompoundSampleFilter::number_of_clauses() const { return m_clauses.size() ; }

	SampleFilter const& CompoundSampleFilter::clause( std::size_t i ) const {
		assert( i < m_clauses.size() ) ;
		return m_clauses[i] ;
	}

	bool SampleFilterDisjunction::test( genfile::CohortIndividualSource const& source, std::size_t sample, DetailBlock* detail ) const {
		bool result = false ;
		for( std::size_t i = 0; i < number_of_clauses(); ++i ) {
			bool a = clause( i ).test( source, sample ) ;
			if( detail ) {
				(*detail)( 0, i ) = a ;
			}
			result = result || a ;
		}
		return result ;
	}

	void SampleFilterDisjunction::summarise( std::ostream& o ) const {
		if( number_of_clauses() > 1 ) {
			o << "( " ;
		}
		for( std::size_t i = 0; i < number_of_clauses(); ++i ) {
			if( i > 0 ) {
				o << " OR " ;
			}
			o << clause( i ) ;
		}
		if( number_of_clauses() > 1 ) {
			o << " )" ;
		}
	}

	bool SampleFilterConjunction::test( genfile::CohortIndividualSource const& source, std::size_t sample, DetailBlock* detail ) const {
		for( std::size_t i = 0; i < number_of_clauses(); ++i ) {
			bool a = clause( i ).test( source, sample ) ;
			if( detail ) {
				(*detail)( 0, i ) = a ;
			}
			if( !a ) {
				return false ;
			}
		}
		return true ;
	}
	
	void SampleFilterConjunction::summarise( std::ostream& o ) const {
		if( number_of_clauses() > 1 ) {
			o << " (" ;
		}
		for( std::size_t i = 0; i < number_of_clauses(); ++i ) {
			if( i > 0 ) {
				o << " AND " ;
			}
			o << clause( i ) ;
		}
		if( number_of_clauses() > 1 ) {
			o << " )" ;
		}
	}
	
}


