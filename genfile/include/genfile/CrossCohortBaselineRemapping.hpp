
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_CROSS_COHORT_BASELINE_REMAPPING_HPP
#define GENFILE_CROSS_COHORT_BASELINE_REMAPPING_HPP

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <boost/bimap.hpp>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/CrossCohortCovariateValueMapping.hpp"

namespace genfile {
	class CrossCohortBaselineRemapping: public CrossCohortCovariateValueMapping {
	public:
		typedef std::auto_ptr< CrossCohortBaselineRemapping > UniquePtr ;
		static UniquePtr create( CrossCohortCovariateValueMapping::UniquePtr source, Entry const baseline ) ;
	public:
		CrossCohortBaselineRemapping( CrossCohortCovariateValueMapping::UniquePtr source, Entry const baseline ) ;
		void add_source( CohortIndividualSource const& source ) ;
		Histogram const& histogram() const ;
		std::size_t get_number_of_unmapped_values() const ;
		std::size_t get_number_of_distinct_unmapped_values() const ;
		std::size_t get_number_of_distinct_mapped_values() const ;
		std::size_t get_number_of_missing_values() const ;
		
		// Get the mapped value of the given entry.
		// Calling this with a value not included among the cohorts passed to this mapping
		// results in undefined behaviour.
		Entry get_mapped_value( Entry const& level ) const ;

		// Get the unmapped value (the value occuring within the cohort sources) corresponding
		// to the given mapped value.  Calling this with a value not included among the cohorts
		// passed to this mapping results in undefined behaviour.
		Entry get_unmapped_value( Entry const& level ) const ;

		// Return a human-readable summary of the mapping and/or its values.
		std::string get_summary( std::string const& prefix = "" ) const ;

	private:
		CrossCohortCovariateValueMapping::UniquePtr m_source ;
		Entry const m_baseline ;
		typedef boost::bimap< Entry, Entry > Map ;
		Map m_map ;
		
	private:
		void recompute_map() ;
	} ;
}

#endif
