#ifndef GENFILE_QUANTILE_NORMALISING_CROSS_COHORT_COVARIATE_MAPPING_HPP
#define GENFILE_QUANTILE_NORMALISING_CROSS_COHORT_COVARIATE_MAPPING_HPP

#include <boost/math/distributions/normal.hpp>
#include "genfile/CrossCohortCovariateValueMapping.hpp"

namespace genfile {
	class QuantileNormalisingCrossCohortCovariateValueMapping: public ContinuousVariableCrossCohortCovariateValueMapping
		// A mapping which scales its given values (which must be real numbers) so that
		// they are approximately normally distributed with mean 0 and variance 1.
	{
	public:
		QuantileNormalisingCrossCohortCovariateValueMapping(
			std::string const& column_name
		) ;
		
		void add_source( CohortIndividualSource const& source ) ;
		std::size_t get_number_of_distinct_mapped_values() const ;
		
		Entry get_unmapped_value( Entry const& level ) const ;
		Entry get_mapped_value( Entry const& entry ) const ;
		
		// std::string get_summary( std::string const& prefix = "" ) const ;
		
		virtual std::string get_mapping_name() const { return "quantile normalised" ; }
		
	private:
		typedef std::map< double, double > Mapping ;
		Mapping m_mapping ;
		Mapping m_reverse_mapping ;
	private:
		void analyse_values( Entries const& entries ) ;
		boost::math::normal m_normal_distribution ;
	} ;
}

#endif
