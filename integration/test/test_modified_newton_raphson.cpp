
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <boost/bind.hpp>
#include <boost/noncopyable.hpp>
#include "test_case.hpp"
#include "integration/ModifiedNewtonRaphson.hpp"
#include <Eigen/Dense>

namespace impl {
	struct FunctionBase: public boost::noncopyable {
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::MatrixXd Matrix ;
	} ;

	struct Function1: public FunctionBase {
		// This function f has four local maxima and one local minimum.
		// We should walk to one of the maxima depending on where we start.
		// For appropriately modified NR, we should never walk to the minimum.
		// f' = -0.1 x^4 - 0.1 y^4 + x^2 + y^2
		// This has a minimum at (0,0)
		// And maxima at nonzero solutions of
		// -0.4 x^3 - 2x = 0 and -0.4 y^3 + 2y = 0
		// (x,y) = (\pm sqrt(5), \pm sqrt(5))
		//
		// First derivative is
		//
		// f' = ( -0.4 x^3 + 2x, -0.4 y^3 + 2y )
		//
		// Diagonals of second derivative are:
		// f'' = -1.2 x^2 + 2, -1.2 y^2 + 2
		//
		// With zeroes on off-diagonals
		//
		void evaluate_at( point, std::size_t number_of_derivatives = 0 ) {
			assert( point.size() == 2 ) ;
			m_point = point ;
		}
		double get_value_of_function() const {
			return -0.1 * ( std::pow( m_point(0), 4 ) + std::pow( m_point(1), 4 ) )
				+ m_point(0) * m_point(0)
				+ m_point(1) * m_point(1)
			;
		}
		Eigen::VectorXd get_value_of_first_derivative() const {
			Eigen::VectorXd result ;
			result.setZero(2) ;
			result(0) = -0.4 * std::pow( m_point(0), 3 ) + 2.0 * m_point(0) ;
			result(1) = -0.4 * std::pow( m_point(1), 3 ) + 2.0 * m_point(1) ;
			return result ;
		}

		Eigen::MatrixXd get_value_of_first_derivative() const {
			Eigen::MatrixXd result ;
			result.setZero(2,2) ;
			result(0,0) = -1.2 * std::pow( m_point(0), 2 ) + 2.0 ;
			result(1,1) = -1.2 * std::pow( m_point(1), 2 ) + 2.0 ;
			return result ;
		}

	private:
		Eigen::Vector m_point ;
	} ;
}

AUTO_TEST_CASE( test_modified_newton_raphson_2d ) {
	{
		impl::Function1 function ;
		
		
	}
}
