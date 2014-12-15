
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <Eigen/Core>
#include "genfile/vcf/get_set_eigen.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"
#include "components/SampleSummaryComponent/IntensityDistributionComputation.hpp"
#include "metro/mean_and_variance.hpp"

// #define DEBUG_INTENSITY_DISTRIBUTION_COMPUTATION 1

namespace sample_stats {
	IntensityDistributionComputation::IntensityDistributionComputation():
		m_call_threshhold( 0.9 ),
		m_snp_index( 0 )
	{}

	void IntensityDistributionComputation::accumulate( genfile::SNPIdentifyingData const&, Genotypes const& genotypes, genfile::VariantDataReader& data_reader ) {
		if( !data_reader.supports( "XY" )) {
			return ;
		}

		int const N = genotypes.rows() ;

		{
			genfile::vcf::MatrixSetter< IntensityMatrix > intensity_setter( m_intensities, m_nonmissingness ) ;
			data_reader.get( "XY", intensity_setter ) ;
			assert( m_intensities.rows() == N ) ;
			assert( m_intensities.cols() == 2 ) ;
			assert( m_nonmissingness.rows() == N ) ;
			assert( m_nonmissingness.cols() == 2 ) ;


			// Explicit check for nan or inf.  These mess with our computations so we exclude 'em.
			// Can't really find a better way to do this than a loop right now.
			for( int i = 0; i < N; ++i ) {
				for( int j = 0; j < 2; ++j ) {
					if( ( m_intensities(i,j) - m_intensities(i,j) ) != ( m_intensities(i,j) - m_intensities(i,j) ) ) {
						m_intensities(i,j) = 0 ;
						m_nonmissingness(i,j) = 0 ;
					}
				}
			}	
		}
		
		if( m_snp_index == 0 ) {
			m_number_of_samples = N ;
			m_accumulator_by_genotype.resize( 4 ) ;
		}

		// We have a column for X and a column for Y.
		// We use a linear transformation to compute means of X, Y, X+Y and X-Y.
		Eigen::Matrix< double, 2, 4 > transform ;
		transform <<
			1,  0,  1,  1,
			0,  1,  1, -1
		;
		m_intensities = m_intensities * transform ;
		m_nonmissingness = m_nonmissingness * transform ;

		// In third column, nonmissingness now equals #of non-missing values (zero, one or two).
		// Replace 3rd and fourth columns with an indicator.
		m_nonmissingness.col(2).array() = ( m_nonmissingness.col(2).array() == 2 ).cast< double >() ;
		m_nonmissingness.col(3) = m_nonmissingness.col(2) ;
		
#if DEBUG_INTENSITY_DISTRIBUTION_COMPUTATION
		std::cerr << "Intensities are:\n"
			<< m_intensities.block( 0, 0, std::min( N, 10 ), 4 )
			<< "with nonmissingness:\n"
			<< m_nonmissingness.block( 0, 0, std::min( N, 10 ), 4 )
			<< "\n" ;
#endif
		// This is an "on-line" algorithm for computing the mean and variance.
		// see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#On-line_algorithm
		//
		
		m_accumulator.accumulate( m_intensities, m_nonmissingness ) ;

#if DEBUG_INTENSITY_DISTRIBUTION_COMPUTATION
			std::cerr << "mean is " << m_accumulator.get_mean().block( 0, 0, std::min( N, 10 ), 2 ) << "\n" ;
#endif
		for( int g = 0; g < 4; ++g ) {
			m_intensities_by_genotype = m_intensities.block( 0, 0, N, 2 ) ;
			m_nonmissingness_by_genotype = m_nonmissingness.block( 0, 0, N, 2 ) ;
			
			for( int i = 0; i < N; ++i ) {
				if( g == 3 ) { // representing missing genotype
					if( genotypes.row( i ).maxCoeff() >= m_call_threshhold ) {
						m_nonmissingness_by_genotype.row( i ).setZero() ;
					}
				}
				else {
					if( genotypes( i, g ) < m_call_threshhold ) {
						m_nonmissingness_by_genotype.row( i ).setZero() ;
					}
				}
			}

#if DEBUG_INTENSITY_DISTRIBUTION_COMPUTATION
			std::cerr << "Accumulating genotypes for g=" << g << ":\n"
				<< m_intensities_by_genotype.block( 0, 0, std::min( N, 10 ), 2 )
				<< "\n nonmissingness = "
				<< m_nonmissingness_by_genotype.block( 0, 0, std::min( N, 10 ), 2 )
				<< "\n" ;
#endif
			m_accumulator_by_genotype[ g ].accumulate( m_intensities_by_genotype, m_nonmissingness_by_genotype ) ;
#if DEBUG_INTENSITY_DISTRIBUTION_COMPUTATION
			std::cerr << "mean is " << m_accumulator_by_genotype[ g ].get_mean().block( 0, 0, std::min( N, 10 ), 2 ) << "\n" ;
#endif
		}

		++m_snp_index ;
	}

	void IntensityDistributionComputation::compute( int sample, ResultCallback callback ) {
		if( m_snp_index > 0 ) {
			Eigen::MatrixXd const& mean = m_accumulator.get_mean() ;
			Eigen::MatrixXd const& variance = m_accumulator.get_variance() ;
			assert( std::size_t( mean.rows() ) == m_number_of_samples ) ;
			assert( std::size_t( mean.cols() ) == 4 ) ;
		
			callback( sample, "mean_X", mean( sample, 0 ) ) ;
			callback( sample, "mean_Y", mean( sample, 1 ) ) ;
			callback( sample, "mean_X+Y", mean( sample, 2 ) ) ;
			callback( sample, "mean_X-Y", mean( sample, 3 ) ) ;
			callback( sample, "variance_X", variance( sample, 0 ) ) ;
			callback( sample, "variance_Y", variance( sample, 1 ) ) ;
			callback( sample, "variance_X+Y", variance( sample, 2 ) ) ;
			callback( sample, "variance_X-Y", variance( sample, 3 ) ) ;

			for( int g = 0; g < 4; ++g ) {
				std::string const stub = "g=" + ( g == 3 ? std::string( "NA" ) : genfile::string_utils::to_string( g ) ) ;
				callback( sample, stub + ":mean_X", m_accumulator_by_genotype[g].get_mean( sample, 0 ) ) ;
				callback( sample, stub + ":mean_Y", m_accumulator_by_genotype[g].get_mean( sample, 1 ) ) ;
				callback( sample, stub + ":variance_X", m_accumulator_by_genotype[g].get_variance( sample, 0 ) ) ;
				callback( sample, stub + ":variance_Y", m_accumulator_by_genotype[g].get_variance( sample, 1 ) ) ;
			}
		}
	}

	std::string IntensityDistributionComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "IntensityDistributionComputation" ;
	}
}

