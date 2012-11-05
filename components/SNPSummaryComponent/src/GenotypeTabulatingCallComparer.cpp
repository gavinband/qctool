
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <map>
#include <cmath>
#include <cassert>
#include <boost/function.hpp>
#include <Eigen/Core>
#include "components/SNPSummaryComponent/GenotypeTabulatingCallComparer.hpp"

GenotypeTabulatingCallComparer::GenotypeTabulatingCallComparer():
	m_threshhold( 0.9 )
{}

void GenotypeTabulatingCallComparer::compare(
	Eigen::MatrixXd const& left,
	Eigen::MatrixXd const& right,
	Callback callback
) const {
	// Make table of counts
	int const N = left.rows() ;
	assert( right.rows() == N ) ;
	assert( left.cols() == 3 ) ;
	assert( right.cols() == 3 ) ;

	Eigen::MatrixXd table( 4, 4 ) ;

	for( std::size_t i = 0; i < N; ++i ) {
		Eigen::MatrixXd::Index left_best_call = -1 ;
		Eigen::MatrixXd::Index right_best_call = -1 ;
		double left_best_call_prob = left.row( i ).maxCoeff( &left_best_call ) ;
		double right_best_call_prob = right.row( i ).maxCoeff( &right_best_call ) ;

		if( left_best_call_prob < m_threshhold ) {
			left_best_call = 3 ;
		}
		if( right_best_call_prob < m_threshhold ) {
			right_best_call = 3 ;
		}
		++table( left_best_call, right_best_call ) ;
	}

	for( int g1 = 0; g1 < 4; ++g1 ) {
		std::string const variable_name1 = ( g1 == 3 ) ? "missing" : genfile::string_utils::to_string( g1 ) ;
		for( int g2 = 0; g2 < 4; ++g2 ) {
			std::string const variable_name2 = ( g2 == 3 ) ? "missing" : genfile::string_utils::to_string( g2 ) ;
			callback( variable_name1 + "/" + variable_name2, table( g1, g2 ) ) ;
		}
	}
}
