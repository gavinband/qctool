
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_NORMALISING_CROSS_COHORT_COVARIATE_MAPPING_HPP
#define GENFILE_NORMALISING_CROSS_COHORT_COVARIATE_MAPPING_HPP

#include <vector>
#include <map>
#include <set>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/ContinuousVariableCrossCohortCovariateValueMapping.hpp"

namespace genfile {
	class NormalisingCrossCohortCovariateValueMapping: public ContinuousVariableCrossCohortCovariateValueMapping
		// A mapping which scales its given values (interpreted as real numbers) so that they have the given
		// mean and variance.
	{
	public:
		NormalisingCrossCohortCovariateValueMapping(
			std::string const& column_name
		) ;
		
		void add_source( CohortIndividualSource const& source ) ;
		
		std::size_t get_number_of_distinct_mapped_values() const ;
		Entry get_unmapped_value( Entry const& level ) const ;
		Entry get_mapped_value( Entry const& entry ) const ;
		
		// std::string get_summary( std::string const& prefix = "" ) const ;
		
		virtual std::string get_mapping_name() const { return "normalised" ; }
		
	private:
		double m_mean, m_variance, m_standard_deviation ;
	private:
		void analyse_values( Histogram const& histogram ) ;
	} ;
}

#endif
