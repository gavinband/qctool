
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
#include "genfile/SampleFilter.hpp"

namespace genfile {
	struct CompoundSampleFilter: public SampleFilter
	{
		typedef std::auto_ptr< CompoundSampleFilter > UniquePtr ;
		
		void add_clause( SampleFilter::UniquePtr clause ) ;
		std::size_t number_of_clauses() const ;
		SampleFilter const& clause( std::size_t i ) const ;
		
	private:
		boost::ptr_vector< SampleFilter > m_clauses ;
	} ;

	struct SampleFilterDisjunction: public CompoundSampleFilter {
		bool test( genfile::CohortIndividualSource const& source, std::size_t sample, DetailBlock* detail ) const ;
		void summarise( std::ostream& ) const ;
	} ;

	struct SampleFilterConjunction: public CompoundSampleFilter {
		bool test( genfile::CohortIndividualSource const& source, std::size_t sample, DetailBlock* detail ) const ;
		void summarise( std::ostream& ) const ;
	} ;
}

#endif

