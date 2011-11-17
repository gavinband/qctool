#include <memory>
#include <Eigen/Core>
#include "Distribution.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "IntensityModel.hpp"
#include "NormalClusterIntensityModel.hpp"

IntensityModel::UniquePtr IntensityModel::estimate( Eigen::MatrixXd const& intensities, genfile::SingleSNPGenotypeProbabilities const& genotypes, double threshhold ) {
	std::size_t number_of_samples = intensities.cols() ;
	assert( genotypes.size() == number_of_samples ) ;
	assert( std::size_t( intensities.cols() ) == number_of_samples ) ;
	assert( intensities.rows() == 2 ) ;

	std::vector< Eigen::Vector2d > means( 3, Eigen::Vector2d::Zero() ) ;
	std::vector< Eigen::Matrix2d > variances( 3, Eigen::Matrix2d::Zero() ) ;
	std::vector< double > non_missing_counts( 3, 0.0 ) ;

	Eigen::Vector2d outlier_mean ;
	double outlier_count = 0.0; 

	for( std::size_t i = 0; i < number_of_samples; ++i ) {
		if( intensities( 0, i ) != intensities( 0, i ) || intensities( 1, i ) != intensities( 1, i ) ) {
			// intensity is missing.  Ignore this point.
		}
		else {
			bool called = false ;
			for( std::size_t g = 0; g < 3; ++g ) {
				if( genotypes( i, g ) >= threshhold ) {
					means[g] += intensities.col( i ) ;
					variances[g] += intensities.col( i ) * intensities.col( i ).transpose() ;
					++non_missing_counts[g] ;
					called = true ;
					break ;
				}
			}
			if( !called ) {
				outlier_mean += intensities.col( i ) ;
				++outlier_count ;
			}
		}
	}

	// complete mean and variance computation.
	for( std::size_t g = 0; g < 3; ++g ) {
		means[g] /= non_missing_counts[g] ;
		variances[g] -= means[g] * means[g].transpose() ;
		variances[g] /= ( non_missing_counts[g] - 1.0 ) ;
	}
	
	outlier_mean /= outlier_count ;
	
	double overall_non_missing_count = outlier_count + std::accumulate( non_missing_counts.begin(), non_missing_counts.end(), 0.0 ) ;
	
	NormalClusterIntensityModel::UniquePtr result( new NormalClusterIntensityModel() ) ;
	for( std::size_t g = 0; g < 3; ++g ) {
		result->add_cluster( means[ g ], variances[ g ], non_missing_counts[g] / overall_non_missing_count ) ;
	}
	
	result->add_cluster( outlier_mean, 5.0 * Eigen::Matrix2d::Identity(), outlier_count / overall_non_missing_count ) ;
	return IntensityModel::UniquePtr( result.release() ) ;
}
