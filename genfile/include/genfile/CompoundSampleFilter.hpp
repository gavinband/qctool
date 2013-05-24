
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef GENFILE_COMPOUND_SAMPLE_FILTER_HPP
#define GENFILE_COMPOUND_SAMPLE_FILTER_HPP

#include <memory>
#include <vector>
#include <set>
#include <boost/ptr_container/ptr_vector.hpp>
#include "genfile/VariantEntry.hpp"
#include "SampleFilter.hpp"

namespace genfile {
	struct CompoundSampleFilter: public SampleFilter
	{
		void add_clause( SampleFilter::UniquePtr clause ) {
			m_clauses.push_back( clause ) ;
		}
		std::size_t number_of_clauses() const { return m_clauses.size() ; }
		SampleFilter const& clause( std::size_t i ) const { assert( i < m_clauses.size() ) ; return m_clauses[i] ; }
	private:
		boost::ptr_vector< SampleFilter > m_clauses ;
	} ;

	struct SampleFilterDisjunction: public CompoundSampleFilter {
		bool test( genfile::CohortIndividualSource const&, std::size_t sample ) const {
			bool result = true ;
			for( std::size_t i = 0; i < number_of_clauses(); ++i ) {
				result = result | clause( i ).test( sample ) ;
			}
			return result ;
		}
	} ;

	struct SampleFilterConjunction: public CompoundSampleFilter {
		bool test( genfile::CohortIndividualSource const&, std::size_t sample ) const {
			for( std::size_t i = 0; i < number_of_clauses(); ++i ) {
				if( !clause( i ).test( sample ) ) {
					return false ;
				}
			}
			return true ;
		}
	} ;
}

#endif

