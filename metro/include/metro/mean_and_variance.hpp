
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_MEAN_AND_VARIANCE
#define METRO_MEAN_AND_VARIANCE

#include <limits>
#include <Eigen/Dense>

namespace metro {
	// Compute mean and variance for a vector.
	// All values must be non-missing (i.e. not NaN.)
	template< typename Data >
	std::pair< double, double > compute_mean_and_variance( Data const& data ) {
		double const mean = data.sum() / data.size() ;
		double variance = std::numeric_limits< double >::quiet_NaN() ;
		if( data.size() > 1 ) {
			variance = ( data.array() - mean ).square().sum() / ( data.size() - 1  ) ;
		}
		return std::make_pair( mean, variance ) ;
	}

	// Compute mean and variance for a vector ignoring missing values.
	template< typename Data >
	std::pair< double, double > compute_mean_and_variance( Data const& data, Data const& nonmissingness ) {
		assert( data.size() == nonmissingness.size() ) ;
		// Ensure all non-missing values are zero.
		double const mean = ( data.array() * nonmissingness.array() ).sum() / nonmissingness.sum() ;
		double variance = std::numeric_limits< double >::quiet_NaN() ;
		if( data.size() > 1 ) {
			variance = (( data.array() - mean ) * nonmissingness.array() ).square().sum() / ( nonmissingness.sum() - 1 ) ;
		}
		return std::make_pair( mean, variance ) ;
	}
	
	
	
	// This struct implements an "on-line" algorithm for computing the mean and variance.
	// see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#On-line_algorithm
	// This implementation computes mean and per-element variance (but not covariance) of a matrix
	// of values.
	struct OnlineElementwiseMeanAndVariance {
		typedef Eigen::MatrixXd Storage ;
	public:
		template< typename Data, typename Nonmissingness >
		void accumulate( Data const& data, Nonmissingness const& nonmissingness ) {
			assert( data.rows() == nonmissingness.rows() ) ;
			assert( data.cols() == nonmissingness.cols() ) ;
			if( m_mean.rows() == 0 ) {
				m_nonmissingness = nonmissingness ;
				// resize storage to match data
				m_mean.setZero( data.rows(), data.cols() ) ;
				m_sum_of_squares_of_differences.setZero( data.rows(), data.cols() ) ;
			} else {
				assert( data.rows() == m_mean.rows() ) ;
				assert( data.cols() == m_mean.cols() ) ;
				assert( nonmissingness.rows() == m_mean.rows() ) ;
				assert( nonmissingness.cols() == m_mean.cols() ) ;
				m_nonmissingness += nonmissingness ;
			}
			m_delta = data - m_mean ;
			
			//std::cerr << "m_nonmissingness =\n" << m_nonmissingness.block( 0, 0, 10, 4 ) << "\n" ;
			//std::cerr << "m_delta =\n" << m_delta.block( 0, 0, 10, 4 ) << "\n" ;
			//std::cerr << "m_mean =\n" << m_mean.block( 0, 0, 10, 4 ) << "\n" ;
			m_mean.array() += nonmissingness.array() * ( m_delta.array() / ( m_nonmissingness.array() + ( m_nonmissingness.array() == 0 ).cast< double >() )) ;
			//std::cerr << "m_mean after update =\n" << m_mean.block( 0, 0, 10, 4 ) << "\n" ;
			m_sum_of_squares_of_differences.array() += nonmissingness.array() * ( m_delta.array() * ( data - m_mean ).array() ) ;
		}
		
		template< typename Data >
		void accumulate( Data const& data ) {
			Storage const nonmissingness = Storage::Constant( data.rows(), data.cols(), 1 ) ;
			this->accumulate( data, nonmissingness ) ;
		}

		Storage get_mean() const ;
		Storage get_variance() const ;
		
	private:
		Storage m_nonmissingness ;
		Storage m_mean ;
		Storage m_sum_of_squares_of_differences ;
		Storage m_delta ;
	} ;
}

#endif
