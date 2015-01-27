//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CORE_MODIFIED_NEWTON_RAPHSON_HPP
#define CORE_MODIFIED_NEWTON_RAPHSON_HPP

#define DEBUG_NEWTON_RAPHSON 1

#include <limits>
#include <iostream>
#include <iomanip>
#include "integration/Error.hpp"

namespace integration {
	namespace {
		template< typename Function, typename DirectionPicker, typename StoppingCondition >
		struct Stepper: public boost::noncopyable {
			Stepper(
				DirectionPicker& solver,
				StoppingCondition& stopping_condition
			):
				 m_solver( solver ),
				 m_stopping_condition( stopping_condition )
			{
			}

			bool step( Function& function, typename Function::Vector const& point, typename Function::Vector* result ) {
				*result = m_solver.compute( function, point ) ;
				return !m_stopping_condition( function, *result ) ;
			}

		private:
			DirectionPicker& m_solver ;
			StoppingCondition& m_stopping_condition ;
		} ;
		
		template< typename Vector >
		double oneNorm( Vector const& v ) {
			return v.array().abs().maxCoeff() ;
		}
	}
	//
	// The following function represents a generic version of the multidimensional modified Newton-Raphson
	// algorithm which finds a local maximum of the given function.
	// It takes an 'function' object which evaluates the function and its first and second derivatives.
	// A 'direction_picker' object is used to pick an ascent direction at each step.
	// If a non-ascent direction is picked, the algorithm steps along it unconditionally without line search.
	// The iteration stops when the picked direction is zero.
	//
	template< typename Function, typename DirectionPicker, typename StoppingCondition >
	typename Function::Vector find_maximum_by_modified_newton_raphson_with_line_search(
		Function& function,
		typename Function::Vector currentPoint,
		DirectionPicker& direction_picker,
		StoppingCondition& stopping_condition,
		double const line_search_tolerance = 0.1
	)
	{
		typedef typename Function::Vector Vector ;
		typedef typename Function::Matrix Matrix ;
		
		// 
		// If h is an ascent direction (f(x + lambda h) > f(x) for small enough lambda)
		// then we can do a line search starting with lambda = 1.
		// We'll have f(x + lambda h) = f(x) + lambda Df_x (h) + O(lambda^2).
		// Can look for a lambda with f(x + lambda h) - f(x) > c * lambda Df_x (h), say,
		// where lambda = the line search tolerance.
		//

		Vector first_derivative ;
		Vector h ;
		Stepper< Function, DirectionPicker, StoppingCondition > stepper( direction_picker, stopping_condition ) ;
		for(
			function.evaluate_at( currentPoint, 1 ) ;
			stepper.step( function, currentPoint, &h ) ;
			function.evaluate_at( currentPoint, 1 )
		) {
			double const current_function_value = function.get_value_of_function() ;
			first_derivative = function.get_value_of_first_derivative() ;
			double directional_derivative = first_derivative.transpose() * h ;

#if DEBUG_NEWTON_RAPHSON
			std::cerr << "MNR: --------------------------\n" ;
			std::cerr << "MNR: point is                 : " << currentPoint.transpose() << ".\n" ;
			std::cerr << "MNR: first derivative is      : " << first_derivative.transpose() << " with 1-norm " << oneNorm( first_derivative ) << ".\n" ;
			std::cerr << "MNR: line search direction is : " << h.transpose() << ".\n" ;
			std::cerr << "MNR: current function value is: " << std::setprecision(25) << current_function_value << ".\n" ;
			std::cerr << "MNR: directional derivative is: " << directional_derivative << ".\n" ;
			std::cerr << "MNR: machine epsilon is: " << std::numeric_limits< double >::epsilon() << ".\n" ;
			std::cerr << "MNR: starting line search...\n" ;
#endif

			double lambda = 1 ;
			// We can't reasonably expect to improve the function if the directional derivative is close to epsilon.
			// But we might still make the derivative smaller.
			if( directional_derivative > ( 10 * std::numeric_limits< double >::epsilon() ) ) {
				// Directional derivative is not too small, so we can hope to improve the function value.
				for(
					function.evaluate_at( currentPoint + lambda * h, 0 ) ;
					( function.get_value_of_function() - current_function_value ) < ( line_search_tolerance * lambda * directional_derivative ) ;
					function.evaluate_at( currentPoint + lambda * h, 0 )
				) {
	#if DEBUG_NEWTON_RAPHSON
				std::cerr << "MNR: tried line search with lambda=" << lambda << ": function = " << std::setprecision(25) << function.get_value_of_function() << ".\n" ;
	#endif
					lambda = lambda * 0.9 ;
					if( lambda < std::numeric_limits< double >::epsilon() ) {
						throw NumericalError(
							"integration::find_maximum_by_modified_newton_raphson_with_line_search()",
							"Line search failed."
						) ;
					}
				}
			} else {
				// Directional derivative is tiny.
				// We can't reasonably expect to improve the function if the directional derivative is close to epsilon.
				// We aim to not make it worse, and to make the derivative smaller.
				double const derivativeOneNorm = oneNorm( first_derivative ) ;
				for(
					function.evaluate_at( currentPoint + lambda * h, 1 ) ;
					(( function.get_value_of_function() - current_function_value ) < 0 )
					|| ( oneNorm( function.get_value_of_first_derivative() ) > derivativeOneNorm ) ;
					function.evaluate_at( currentPoint + lambda * h, 0 )
				) {
	#if DEBUG_NEWTON_RAPHSON
				std::cerr << "MNR: tried line search with lambda=" << lambda
					<< ": function = " << std::setprecision(25) << function.get_value_of_function()
					<< ": derivative = " << function.get_value_of_first_derivative().transpose()
					<< ": 1-norm = " << oneNorm( function.get_value_of_first_derivative() )
					<< ".\n" ;
	#endif
					lambda = lambda * 0.5 ;
					if( lambda < std::numeric_limits< double >::epsilon() ) {
						throw NumericalError(
							"integration::find_maximum_by_modified_newton_raphson_with_line_search()",
							"Line search failed."
						) ;
					}
				}
			}
#if DEBUG_NEWTON_RAPHSON
		std::cerr << "MNR: line search complete with lambda = "
				<< lambda
				<< " and function value "
				<< function.get_value_of_function()
				<< ", improvement = "
				<< ( function.get_value_of_function() - current_function_value )
				<< ".\n" ;
#endif
			
			currentPoint += lambda * h ;
		}
#if DEBUG_NEWTON_RAPHSON
		std::cerr << "MNR: point = " << currentPoint << ", function = " << function.get_value_of_function() << ", first derivative = " << function.get_value_of_first_derivative() << ".\n" ;
#endif
		return currentPoint ;
	}
}

#endif
