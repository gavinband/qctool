//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_MODIFIED_NEWTON_RAPHSON_HPP
#define METRO_MODIFIED_NEWTON_RAPHSON_HPP

#include <limits>
#include <iostream>
#include <iomanip>
#include "Error.hpp"

// #define DEBUG_NEWTON_RAPHSON 1

namespace metro {
	namespace {
		// compute the one-norm
		template< typename Vector >
		double oneNorm( Vector const& v ) {
			return v.array().abs().maxCoeff() ;
		}
	}

	//
	// The following function implements the
	// multidimensional modified Newton-Raphson algorithm.
	// It finds a local maximum of the given function.
	// It takes a 'function' object which evaluates the function and its
	// first and second derivatives.

	// A 'stepper' object is used to tell the algorithm whether to take another
	// step, and if so in which direction.	This should be callable as
	//	`stepper.step( function, current_point, &search_direction )`
	// A false return value indicates the algorithm should stop.	Otherwise, the
	// function should population the search direction and return true.
	
	// If an ascent direction (=directional derivative positive) is picked, the
	// algorithm performs a backward line search to find a step that actually
	// increases the function value.
	// If with nonpositive directional derivative is picked then the algorithm
	// gamely tries to improve the function or derivative anyway, but will bail
	// out and return the current value
	// if it cannot make any improvement.
	//
	template< typename Function, typename Stepper >
	typename Function::Vector find_maximum_by_modified_newton_raphson_with_line_search(
		Function& function,
		typename Function::Vector currentPoint,
		Stepper& stepper,
		double const line_search_tolerance = 0.1
	)
	{
		typedef typename Function::Vector Vector ;
		typedef typename Function::Matrix Matrix ;
		double const epsilon = std::numeric_limits< double >::epsilon() ;
		
		Vector first_derivative ;
		Vector h ;
		for(
			function.evaluate_at( currentPoint, 1 ) ;
//			function.evaluate_at( currentPoint, 1 ) ; // needs likelihood optimisation
			stepper.step( function, currentPoint, &h ) ;
			function.evaluate_at( currentPoint, 1 )
//			function.evaluate_at( currentPoint, 1 ) // needs likelihood optimisation
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
			std::cerr << "MNR: current       || df || is: " << std::setprecision(25) << (first_derivative.array().abs().maxCoeff()) << ".\n" ;
			std::cerr << "MNR: directional derivative is: " << directional_derivative << ".\n" ;
			std::cerr << "MNR: epsilon is: " << epsilon << ".\n" ;
			std::cerr << "MNR: starting line search...\n" ;
#endif
			// Conduct a line search, starting at x + lh with l = 1.
			// We stop if either the Armijo condition is satisfied, i.e. stopping the search when
			//   f(x+lh) > c * l Df_x(h)
			// where c is the line search tolerance.
			
			// OR if the directional derivative is small, we stop if the 
			// function does not decrease more than a small multiple of epsilon
			// and the derivative decreases.
			double lambda = 1 ;
			// We reckon there's no point trying to use the armijo condition if the directional derivative is smaller than this scale
			bool const useArmijo = (( directional_derivative / std::abs( current_function_value )) > 10 * epsilon ) ;
			bool const useDerivative = !useArmijo ;
			double const derivativeOneNorm = oneNorm( first_derivative ) ;
			for(
				function.evaluate_at( currentPoint + lambda * h, useDerivative ? 1 : 0 ) ;
				!(
					(
						useArmijo
						&& ( function.get_value_of_function() - current_function_value ) > ( line_search_tolerance * lambda * directional_derivative ))
					|| 
					(
						useDerivative
						&& (oneNorm( function.get_value_of_first_derivative() ) < derivativeOneNorm)
						&& (( function.get_value_of_function() - current_function_value ) / std::abs( current_function_value )) > -(epsilon*1E6)
					)
				);
				function.evaluate_at( currentPoint + lambda * h, useDerivative ? 1 : 0 )
			) {
#if DEBUG_NEWTON_RAPHSON
			std::cerr << "MNR: tried line search (" << ( useArmijo ? "armijo" : "derivative" ) << ") with lambda=" << lambda
				<< ":      f = " << std::setprecision(25) << function.get_value_of_function()
				<< ";   diff = " << std::setprecision(25) << (function.get_value_of_function() - current_function_value) << "\n"
				<< "; ||df|| = " << std::setprecision(25) << (function.get_value_of_first_derivative().array().abs().maxCoeff())
				<< ".\n" ;
#endif
				lambda *= 0.875 ; // (1-1/8), exactly representable in float
				if( (function.get_value_of_function() - current_function_value) == 0 || lambda < 0.1 ) {
					lambda = 0 ;
					break ;
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
		std::cerr << "MNR: point = " << currentPoint.transpose() << ", function = " << function.get_value_of_function()
			<< ",\nfirst derivative = " << function.get_value_of_first_derivative().transpose() << ".\n" ;
#endif
		return currentPoint ;
	}
}

#endif
