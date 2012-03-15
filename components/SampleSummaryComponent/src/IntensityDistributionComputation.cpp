#ifndef QCTOOL_INTENSITY_DISTRIBUTION_COMPUTATION_HPP
#define QCTOOL_INTENSITY_DISTRIBUTION_COMPUTATION_HPP

#include <Eigen/Core>
#include "genfile/vcf/get_set_eigen.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"

namespace sample_stats {
	IntensityDistributionComputation::IntensityDistributionComputation():
		m_snp_index( 0 )
	{}

	void IntensityDistributionComputation::accumulate( genfile::SNPIdentifyingData const&, Genotypes const&, genfile::VariantDataReader& ) {
		{
			genfile::vcf::MatrixSetter< IntensityMatrix > intensity_setter( m_intensities ) ;
			data_reader.get( "XY", intensity_setter ) ;
		}
		
		if( m_snp_index++ == 0 ) {
			m_means.setZero( 2, m_intensities.cols() ) ;
			m_variances.setZero( 2, m_intensities.cols() ) ;
			m_sum_of_squares_of_differences.setZero( 2, m_intensities.cols() ) ;
		}

		// see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#On-line_algorithm
		m_difference = intensities - m_means ;
		m_means += m_difference / m_snp_index ;
		m_sum_of_squares_of_differences += m_difference.array() * ( intensities - m_means ).array() ;
	}

	void IntensityDistributionComputation::compute( ResultCallback ) {
		m_variances = m_sum_of_squares_of_differences / ( m_snp_index + 1 ) ;
	}

	std::string IntensityDistributionComputation::get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
		
	}
}

#endif
