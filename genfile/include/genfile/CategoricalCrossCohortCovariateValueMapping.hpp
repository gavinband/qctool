
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_CATEGORICAL_CROSS_COHORT_COVARIATE_MAPPING_HPP
#define GENFILE_CATEGORICAL_CROSS_COHORT_COVARIATE_MAPPING_HPP

#include <vector>
#include <map>
#include <set>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/LevelCountingCrossCohortCovariateValueMapping.hpp"

namespace genfile {
	class CategoricalCrossCohortCovariateValueMapping: public LevelCountingCrossCohortCovariateValueMapping
		// A mapping which maps its values onto the inclusive range of positive integers
		// 1...number of distinct values.
	{
	public:
		CategoricalCrossCohortCovariateValueMapping(
			std::string const& column_name
		) ;
		
		std::size_t get_number_of_distinct_mapped_values() const ;

		// Return the entry corresponding to level i.
		// i must be a positive integer in the range 1...number of distinct entries
		Entry get_unmapped_value( Entry const& level ) const ;
		
		// Return the level (positive integer) corresponding to the given entry.
		// It is an error to pass this method an entry not present in the relevant column of the sources
		// this object has had added.
		Entry get_mapped_value( Entry const& entry ) const ;
		
		std::string get_summary( std::string const& prefix = "" ) const ;
		
	} ;
}

#endif
