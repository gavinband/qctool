
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <Eigen/Core>
#include "genfile/vcf/get_set_eigen.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"
#include "components/SampleSummaryComponent/IntensityDistributionComputation.hpp"

namespace sample_stats {
	IntensityDistributionComputation::IntensityDistributionComputation():
		m_snp_index( 0 )
	{}

	void IntensityDistributionComputation::accumulate( genfile::SNPIdentifyingData const&, Genotypes const& genotypes, genfile::VariantDataReader& data_reader ) {
		if( !data_reader.supports( "XY" )) {
			return ;
		}

		{
			genfile::vcf::MatrixSetter< IntensityMatrix > intensity_setter( m_intensities ) ;
			data_reader.get( "XY", intensity_setter ) ;
		}
		
		if( m_snp_index == 0 ) {
			m_means.setZero( 4, m_intensities.cols() ) ;
			m_means2.setZero( 4, m_intensities.cols() ) ;
			m_variances.setZero( 4, m_intensities.cols() ) ;
			m_sum_of_squares_of_differences.setZero( 4, m_intensities.cols() ) ;
		}

		// We have a row for X and a row for Y.
		// Convert to X+Y and X-Y rows.
		Eigen::Matrix< double, 4, 2 > transform ;
		transform <<
			1,  0,
			0,  1,
			1,  1,
			1, -1
		;
		m_intensities = transform * m_intensities ;
		
		// This is an "on-line" algorithm for computing the mean and variance.
		// see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#On-line_algorithm
		//
		++m_snp_index ;
		m_difference = m_intensities - m_means ;
		m_means += m_difference / m_snp_index ;
		m_sum_of_squares_of_differences.array() += m_difference.array() * ( m_intensities - m_means ).array() ;
	}

	void IntensityDistributionComputation::compute( ResultCallback callback ) {
		if( m_snp_index > 0 ) {
			m_variances = m_sum_of_squares_of_differences / ( m_snp_index - 1 ) ;
			for( int sample = 0; sample < m_means.cols(); ++sample ) {
				callback( sample, "X mean", m_means( 0, sample ) ) ;
				callback( sample, "Y mean", m_means( 1, sample ) ) ;
				callback( sample, "X+Y mean", m_means( 2, sample ) ) ;
				callback( sample, "X-Y mean", m_means( 3, sample ) ) ;
				callback( sample, "X variance", m_variances( 0, sample ) / ( m_snp_index - 1 )) ;
				callback( sample, "Y variance", m_variances( 1, sample ) / ( m_snp_index - 1 )) ;
				callback( sample, "X+Y variance", m_variances( 2, sample ) / ( m_snp_index - 1 )) ;
				callback( sample, "X-Y variance", m_variances( 3, sample ) / ( m_snp_index - 1 )) ;
			}
		}
	}

	std::string IntensityDistributionComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "IntensityDistributionComputation" ;
	}
}

