
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_MEAN_AND_COVARIANCE
#define METRO_MEAN_AND_COVARIANCE

#include <limits>

namespace metro {
	// compute mean and covariance of the columns of data.
	// missing data is ignored for this purpose.
	template< typename Data, typename RowVector, typename CovarianceMatrix >
	void compute_mean_and_covariance( Data const& data, Data const& nonmissingness, RowVector& mean, CovarianceMatrix& covariance ) {
		assert( data.cols() == nonmissingness.cols() ) ;
		assert( data.rows() == nonmissingness.rows() ) ;
		typedef typename Data::Index Index ;
		Index const N = data.rows() ;

		// Compute mean
		mean = (( data.array() * nonmissingness.array() ).colwise().sum().array() ) / ( nonmissingness.colwise().sum().array() ) ;

		// compute covariance
		covariance.resize( data.cols(), data.cols() ) ;
		if( data.rows() > 1 ) {
			covariance.setZero() ;
			CovarianceMatrix pairwise_non_missing_counts( data.cols(), data.cols() ) ;
			pairwise_non_missing_counts.setZero() ;
			RowVector mean_centred_row ;
			for( Index i = 0; i < data.rows(); ++i ) {
				mean_centred_row = (( data.row(i) - mean ).array() * nonmissingness.row( i ).array() ) ;
				covariance += mean_centred_row.transpose() * mean_centred_row ;
				pairwise_non_missing_counts += nonmissingness.row( i ).transpose() * nonmissingness.row( i ) ;
			}
			
			// Set any counts less than 2 to NaN.
			pairwise_non_missing_counts.array() *= pairwise_non_missing_counts.array() / pairwise_non_missing_counts.array() ;
			pairwise_non_missing_counts.array() *= ( pairwise_non_missing_counts.array() - 1 ) / ( pairwise_non_missing_counts.array() - 1 ) ;

			covariance.array() /= ( pairwise_non_missing_counts.array() - 1 ) ;
		}
		else {
			covariance.setConstant( std::numeric_limits< double >::quiet_NaN() ) ;
		}
	}

	// Compute mean and covariance, assuming no missing values.
	template< typename Data, typename RowVector, typename CovarianceMatrix >
	void compute_mean_and_covariance( Data const& data, RowVector& mean, CovarianceMatrix& covariance ) {
		typedef typename Data::Index Index ;
		Index const N = data.rows() ;
		// Compute mean
		mean = data.colwise().sum() / N ;
		// compute covariance
		covariance.resize( data.cols(), data.cols() ) ;
		if( data.rows() > 1 ) {
			covariance.setZero() ;
			RowVector mean_centred_row ;
			for( Index i = 0; i < data.rows(); ++i ) {
				mean_centred_row = ( data.row(i) - mean ).array() ;
				covariance += mean_centred_row.transpose() * mean_centred_row ;
			}
			covariance /= ( data.rows() - 1 ) ;
		}
		else {
			covariance.setConstant( std::numeric_limits< double >::quiet_NaN() ) ;
		}
	}
}

#endif
