
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <math.h>
#include <cassert>
#include <Eigen/Core>
#include "metro/log_sum_exp.hpp"

namespace metro {
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
	