
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_CROSS_COHORT_COVARIATE_MAPPING_HPP
#define GENFILE_CROSS_COHORT_COVARIATE_MAPPING_HPP

#include <vector>
#include <map>
#include <set>
#include "genfile/CohortIndividualSource.hpp"

namespace genfile {
	class CrossCohortCovariateValueMapping
		// Base classes for classes which map a given column of phenotypes or covariates,
		// across several cohorts, to transformed values.
	{
	public:
		typedef std::auto_ptr< CrossCohortCovariateValueMapping > UniquePtr ;
		static UniquePtr create( CohortIndividualSource::SingleColumnSpec const& column_spec ) ;

		typedef CohortIndividualSource::Entry Entry ;
		typedef std::map< Entry, unsigned int > Histogram ;
	public:
		virtual ~CrossCohortCovariateValueMapping() {}
		virtual void add_source( CohortIndividualSource const& source ) = 0 ;
		virtual Histogram const& histogram() const = 0 ;
		virtual std::size_t get_number_of_unmapped_values() const = 0 ;
		virtual std::size_t get_number_of_distinct_unmapped_values() const = 0 ;
		virtual std::size_t get_number_of_distinct_mapped_values() const = 0 ;
		virtual std::size_t get_number_of_missing_values() const = 0 ;
		
		// Get the mapped value of the given entry.
		// Calling this with a value not included among the cohorts passed to this mapping
		// results in undefined behaviour.
		virtual Entry get_mapped_value( Entry const& level ) const = 0 ;

		// Get the unmapped value (the value occuring within the cohort sources) corresponding
		// to the given mapped value.  Calling this with a value not included among the cohorts
		// passed to this mapping results in undefined behaviour.
		virtual Entry get_unmapped_value( Entry const& level ) const = 0 ;

		// Return a human-readable summary of the mapping and/or its values.
		virtual std::string get_summary( std::string const& prefix = "" ) const = 0 ;
	} ;
}

#endif
