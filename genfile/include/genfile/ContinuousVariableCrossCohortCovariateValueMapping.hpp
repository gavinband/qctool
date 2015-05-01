
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_CONTINUOUS_VARIABLE_CROSS_COHORT_COVARIATE_MAPPING_HPP
#define GENFILE_CONTINUOUS_VARIABLE_CROSS_COHORT_COVARIATE_MAPPING_HPP

#include <vector>
#include <map>
#include <set>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/LevelCountingCrossCohortCovariateValueMapping.hpp"

namespace genfile {
	class ContinuousVariableCrossCohortCovariateValueMapping: public LevelCountingCrossCohortCovariateValueMapping
	{
		// An identity mapping.  This is only useful insofar as it provides a nice summary of its data.
	public:
		ContinuousVariableCrossCohortCovariateValueMapping(
			std::string const& column_name
		) ;
		
		std::size_t get_number_of_distinct_mapped_values() const ;
		Entry get_unmapped_value( Entry const& level ) const ;
		Entry get_mapped_value( Entry const& entry ) const ;
		
		std::string get_summary( std::string const& prefix = "" ) const ;

	public:
		virtual std::string get_mapping_name() const { return "normalised" ; }

	protected:
		
		static void calculate_mean_and_variance( Histogram const& histogram, double* mean, double* variance ) ;
		typedef std::map< std::pair< double, double >, double > NumericalHistogram ;
		NumericalHistogram get_numerical_histogram(
			Histogram const& histogram,
		 	std::size_t number_of_bins
		) const ;
		std::string print_numerical_histogram(
			NumericalHistogram const& histogram,
			std::string const& prefix,
			std::size_t const height,
			std::size_t const bin_width
		) const ;
	} ;
}

#endif
