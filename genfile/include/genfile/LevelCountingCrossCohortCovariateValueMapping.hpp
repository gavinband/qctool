
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_LEVEL_COUNTING_CROSS_COHORT_COVARIATE_MAPPING_HPP
#define GENFILE_LEVEL_COUNTING_CROSS_COHORT_COVARIATE_MAPPING_HPP

#include <vector>
#include <map>
#include <set>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/CrossCohortCovariateValueMapping.hpp"

namespace genfile {
	class LevelCountingCrossCohortCovariateValueMapping: public CrossCohortCovariateValueMapping
	{
	public:
		LevelCountingCrossCohortCovariateValueMapping(
			std::string const& column_name
		) ;
		
		void add_source( CohortIndividualSource const& source ) ;
		Histogram const& histogram() const ;
		std::size_t get_number_of_unmapped_values() const ;
		std::size_t get_number_of_distinct_unmapped_values() const ;
		std::size_t get_number_of_missing_values() const ;
		
	private:
		std::string m_column_name ;
		Histogram m_histogram ;
		std::size_t m_number_of_missing_values ;
		
	private:

		void add_entries_from_source( CohortIndividualSource const& source, std::string const& column_name ) ;
		// forbid copying / assignment
		LevelCountingCrossCohortCovariateValueMapping( LevelCountingCrossCohortCovariateValueMapping const& other ) ;
		LevelCountingCrossCohortCovariateValueMapping& operator=( LevelCountingCrossCohortCovariateValueMapping const& other ) ;
	} ;
}

#endif

