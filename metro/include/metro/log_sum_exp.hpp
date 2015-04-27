
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_LOG_SUM_EXP_HPP
#define METRO_LOG_SUM_EXP_HPP

#include <Eigen/Core>

namespace metro {
	// Compute the log of a sum of exponentials via an algorithm that 
	template< typename Matrix, typename Nonmissingness >
	double log_sum_exp( Matrix const& data, Nonmissingness const& nonmissingness ) {
		assert( data.rows() == nonmissingness.rows() ) ;
		assert( data.cols() == nonmissingness.cols() ) ;
		if( data.size() == 0 || nonmissingness.sum() == 0 ) {
			return 0.0 ;
		}
		double max_value = ( data.array() * nonmissingness.array() ).maxCoeff() ;
		if( max_value == -std::numeric_limits< double >::infinity() ) {
			return max_value ;
		}
		Eigen::MatrixXd const exponential = ( data - Eigen::MatrixXd::Constant( data.rows(), data.cols(), max_value ) ).array().exp() ;
		return max_value + std::log( ( exponential.array() * nonmissingness.array() ).sum() ) ;
	}

	void rowwise_log_sum_exp( Eigen::MatrixXd const&, Eigen::MatrixXd const& nonmissingness, Eigen::VectorXd* result ) ;
}

#endif
	