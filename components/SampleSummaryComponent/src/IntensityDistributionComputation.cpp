#include <Eigen/Core>
#include "genfile/vcf/get_set_eigen.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"
#include "components/SampleSummaryComponent/IntensityDistributionComputation.hpp"

namespace sample_stats {
	IntensityDistributionComputation::IntensityDistributionComputation():
		m_snp_index( 0 )
	{}

	void IntensityDistributionComputation::accumulate( genfile::SNPIdentifyingData const&, Genotypes const&, genfile::VariantDataReader& data_reader ) {
		if( !data_reader.supports( "XY" )) {
			return ;
		}

		{
			genfile::vcf::MatrixSetter< IntensityMatrix > intensity_setter( m_intensities ) ;
			data_reader.get( "XY", intensity_setter ) ;
		}
		
		if( m_snp_index++ == 0 ) {
			m_means.setZero( 2, m_intensities.cols() ) ;
			m_variances.setZero( 2, m_intensities.cols() ) ;
			m_sum_of_squares_of_differences.setZero( 2, m_intensities.cols() ) ;
		}
		
		// We have a row for X and a row for Y.
		// Convert to X+Y and X-Y rows.
		Eigen::Matrix< double, 2, 2 > transform ;
		transform <<
			1,  1,
			1, -1
		;
		m_intensities = transform * m_intensities ;
		
		// This is an "on-line" algorithm for computing the mean and variance.
		// see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#On-line_algorithm
		//
		m_difference = m_intensities - m_means ;
		m_sum_of_squares_of_differences.array() += m_difference.array() * ( m_intensities - m_means ).array() ;
		m_means += m_difference / m_snp_index ;
	}

	void IntensityDistributionComputation::compute( ResultCallback callback ) {
		if( m_snp_index > 0 ) {
			m_variances = m_sum_of_squares_of_differences / ( m_snp_index - 1 ) ;
			for( int sample = 0; sample < m_means.cols(); ++sample ) {
				callback( sample, "X+Y mean", m_means( 0, sample ) ) ;
				callback( sample, "X-Y mean", m_means( 1, sample ) ) ;
				callback( sample, "X+Y variance", m_variances( 0, sample ) ) ;
				callback( sample, "X-Y variance", m_variances( 1, sample ) ) ;
			}
		}
	}

	std::string IntensityDistributionComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "IntensityDistributionComputation" ;
	}
}

