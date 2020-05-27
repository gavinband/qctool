
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <math.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <Eigen/Core>
#include "metro/log_sum_exp.hpp"

namespace metro {
	template<>
	double log_sum_exp( std::vector< double > const& data ) {
		if( data.size() == 0 ) {
			return -std::numeric_limits< double >::infinity() ;
		}
		std::vector< double >::const_iterator max_i = std::max_element( data.begin(), data.end() ) ;
		double const max_value = *max_i ;
		if( max_value == -std::numeric_limits< double >::infinity() ) {
			return max_value ;
		}
		// exponentiate
		double result = max_value ;
		for( std::size_t i = 0; i < data.size(); ++i ) {
			result += std::exp(data[i] - max_value) ;
		}
		return max_value + std::log( result ) ;
	}
	
	void rowwise_log_sum_exp( Eigen::MatrixXd const& data, Eigen::VectorXd* result ) {
		assert( result ) ;
		assert( result->size() == data.rows() ) ;

		for( int i = 0; i < result->size(); ++i ) {
			(*result)(i) = log_sum_exp( data.row(i) ) ;
		}
	}
	void rowwise_log_sum_exp( Eigen::MatrixXd const& data, Eigen::MatrixXd const& nonmissingness, Eigen::VectorXd* result ) {
		assert( result ) ;
		assert( result->size() == data.rows() ) ;
		assert( data.rows() == nonmissingness.rows() ) ;
		assert( data.cols() == nonmissingness.cols() ) ;

		for( int i = 0; i < result->size(); ++i ) {
			(*result)(i) = log_sum_exp( data.row(i), nonmissingness.row(i) ) ;
		}
	}
}
