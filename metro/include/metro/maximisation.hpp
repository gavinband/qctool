
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CORE_MAXIMISATION_HPP
#define CORE_MAXIMISATION_HPP

#include "metro/NewtonRaphson.hpp"

namespace metro {
	namespace impl {
		// This object represents the derivative of a twice-differentiable function.
		template< typename Function >
		struct Derivative
		{
			typedef typename Function::Vector Vector ;
			typedef typename Function::Matrix Matrix ;

			Derivative( Function& function ): m_function( function ) {}
			Derivative( Derivative const& other ): m_function( other.m_function ) {}

			void evaluate_at( Vector const& parameters ) { m_function.evaluate_at( parameters ) ; }
			Vector get_value_of_function() const { return m_function.get_value_of_first_derivative() ; }
			Matrix get_value_of_first_derivative() const { return m_function.get_value_of_second_derivative() ; }
		private:
			Function& m_function ;
			Derivative& operator=( Derivative const& other ) ;
		} ;

		template< typename Function >
		Derivative< Function > derivative( Function& function ) {
			return Derivative< Function >( function ) ;
		}
	}

	template< typename SmoothFunction, typename Domain, typename StoppingCondition, typename MatrixSolver >
	Domain maximise_by_newton_raphson(
		SmoothFunction& f,
		Domain const& initial_point,
		StoppingCondition& stopping_condition,
		MatrixSolver& solver
	) {
		impl::Derivative< SmoothFunction > Df = impl::derivative( f ) ;
		return find_root_by_newton_raphson( Df, initial_point, stopping_condition, solver ) ;
	}
}

#endif
