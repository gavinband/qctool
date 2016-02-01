
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_INTENSITY_DISTRIBUTION_COMPUTATION_HPP
#define QCTOOL_INTENSITY_DISTRIBUTION_COMPUTATION_HPP

#include <Eigen/Core>
#include "metro/mean_and_variance.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"

namespace sample_stats {
	struct IntensityDistributionComputation: public SampleSummaryComputation
	{
		IntensityDistributionComputation() ;
		void accumulate( genfile::VariantIdentifyingData const&, Genotypes const&, genfile::VariantDataReader& ) ;
		void compute( int sample, ResultCallback ) ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
	private:
		double const m_call_threshhold ;
		typedef Eigen::MatrixXd IntensityMatrix ;
		IntensityMatrix m_intensities ;
		IntensityMatrix m_nonmissingness ;
		IntensityMatrix m_intensities_by_genotype ;
		IntensityMatrix m_nonmissingness_by_genotype ;

		std::size_t m_number_of_samples ;
		std::size_t m_snp_index ;
		metro::OnlineElementwiseMeanAndVariance m_accumulator ;
		std::vector< metro::OnlineElementwiseMeanAndVariance > m_accumulator_by_genotype ;
	} ;
}

#endif
