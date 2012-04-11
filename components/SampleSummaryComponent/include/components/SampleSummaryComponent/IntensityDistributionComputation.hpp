
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_INTENSITY_DISTRIBUTION_COMPUTATION_HPP
#define QCTOOL_INTENSITY_DISTRIBUTION_COMPUTATION_HPP

#include <Eigen/Core>
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"

namespace sample_stats {
	struct IntensityDistributionComputation: public SampleSummaryComputation
	{
		IntensityDistributionComputation() ;
		void accumulate( genfile::SNPIdentifyingData const&, Genotypes const&, genfile::VariantDataReader& ) ;
		void compute( ResultCallback ) ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
	private:
		typedef Eigen::MatrixXd IntensityMatrix ;
		IntensityMatrix m_intensities ;

		std::size_t m_snp_index ;
		IntensityMatrix m_means ;
		IntensityMatrix m_difference ;
		IntensityMatrix m_sum_of_squares_of_differences ;
		IntensityMatrix m_variances ;
	} ;
}

#endif
