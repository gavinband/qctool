//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CORE_MODIFIED_NEWTON_RAPHSON_HPP
#define CORE_MODIFIED_NEWTON_RAPHSON_HPP

#define DEBUG_NEWTON_RAPHSON 1

#include <iostream>
#include <iomanip>
#include "integration/Error.hpp"

namespace integration {
	//
	// The following function represents a generic version of the multidimensional modified Newton-Raphson
	// algorithm which finds a local maximum of the given function.
	// It takes an 'evaluator' object which evaluates the function and its first and second derivatives.
	// A 'solver' object is used to pick a (hopefully ascent) direction at each step.
	// The 'stopping_condition' argument is used to determine when to stop iterating.
	//
	template< typename FunctionEvaluator, typename StoppingCondition, typename MatrixSolver >
	typename FunctionEvaluator::Vector find_maximum_by_modified_newton_raphson_with_line_search(
		FunctionEvaluator& evaluator,
		typename FunctionEvaluator::Vector point,
		StoppingCondition& stopping_condition,
		MatrixSolver& solver,
		double const line_search_tolerance = 0.1
	)
	{
		typedef typename FunctionEvaluator::Vector Vector ;
		typedef typename FunctionEvaluator::Matrix Matrix ;
		evaluator.evaluate_at( point ) ;
		Matrix second_derivative ;
		double function_value = evaluator.get_value_of_function() ;
		Vector first_derivative = evaluator.get_value_of_first_derivative() ;
		Vector h ;

#if DEBUG_NEWTON_RAPHSON
		std::cerr << "MNR: point = " << point << ", function = " << function_value << ", first derivative = " << first_derivative << ".\n" ;
#endif
		while( !stopping_condition( point, function_value, first_derivative ) ) {
			// The Newton-Raphson rule comes from the observation that if
			// f( x + h ) = f( x ) + (D_x f) (h) + higher order terms
			// and if f( x + h ) = 0
			// then h must satisfy (D_x f) (h) = -f( x ) + higher order terms.
			// At each step we solve this and move to the point x + h.
			// If the function is linear, this will actually get us to the root.
			// For modified Newton-Raphson, we instead pick a direction by some
			// modified strategy based on the second (encapsulated in the solver argument).
			// 
			// If h is an ascent direction (f(x + lambda h) > f(x) for small enough lambda)
			// then we can do a line search starting with lambda = 1.
			// We'll have f(x + lambda h) = f(x) + lambda Df_x (h) + O(lambda^2).
			// Can look for a lambda with f(x + lambda h) - f(x) > c * lambda Df_x (h), say,
			// where lambda = the line search tolerance.
			//

#if DEBUG_NEWTON_RAPHSON
			std::cerr << "MNR: --------------------------\n" ;
			std::cerr << "MNR: point is                 : " << point.transpose() << ".\n" ;
			std::cerr << "MNR: first derivative is      : " << first_derivative.transpose() << ".\n" ;
#endif
			second_derivative = evaluator.get_value_of_second_derivative() ;
			solver.compute( second_derivative ) ;
			h = solver.solve( -first_derivative ) ;

			double const current_function_value = evaluator.get_value_of_function() ;
			double lambda = 1 ;
			double directional_derivative = first_derivative.transpose() * h ;

#if DEBUG_NEWTON_RAPHSON
			std::cerr << "MNR: line search direction is : " << h.transpose() << ".\n" ;
			std::cerr << "MNR: current function value is: " << std::setprecision(15) << current_function_value << ".\n" ;
			std::cerr << "MNR: directional derivative is: " << directional_derivative << ".\n" ;
			std::cerr << "MNR: starting line search...\n" ;
#endif
			for(
				evaluator.evaluate_at( point + lambda * h, 0 ) ;
				( evaluator.get_value_of_function() - current_function_value ) < ( line_search_tolerance * lambda * directional_derivative ) ;
				evaluator.evaluate_at( point + lambda * h, 0 )
			) {
#if DEBUG_NEWTON_RAPHSON
				std::cerr << "MNR: tried line search with lambda=" << lambda << ": function = " << evaluator.get_value_of_function() << ".\n" ;
				std::cerr << "MNR: difference = " << ( evaluator.get_value_of_function() - current_function_value ) << ", target = " << ( line_search_tolerance * lambda * directional_derivative ) << ".\n" ;
#endif
				lambda = lambda * 0.5 ;
				if( lambda < 0.0000000001 ) {
					throw NumericalError(
						"integration::find_maximum_by_modified_newton_raphson_with_line_search()",
						"Line search failed."
					) ;
				}
			}
			
#if DEBUG_NEWTON_RAPHSON
			std::cerr << "MNR: line search complete with lambda = "
					<< lambda
					<< " and function value "
					<< evaluator.get_value_of_function()
					<< ", improvement = "
					<< ( evaluator.get_value_of_function() - current_function_value )
					<< ".\n" ;
#endif
			point += lambda * h ;
			evaluator.evaluate_at( point ) ;
#if DEBUG_NEWTON_RAPHSON
			std::cerr << "MNR: point = " << point.transpose() << ", function = " << function_value << ", first derivative = " << first_derivative.transpose() << ".\n" ;
#endif
			first_derivative = evaluator.get_value_of_first_derivative() ;
		}
		return point ;
	}
}

#endif
