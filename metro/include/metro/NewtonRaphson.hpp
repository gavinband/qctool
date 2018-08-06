
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CORE_NEWTON_RAPHSON_HPP
#define CORE_NEWTON_RAPHSON_HPP

// #define DEBUG_NEWTON_RAPHSON 1

namespace metro {
	//
	// The following function represents a generic version of the multidimensional Newton-Raphson
	// algorithm.  It takes an "evaluator" object whose job is to evaluate the function
	// and its derivative simultaneously.  This scheme is useful for functions (such as loglikelihoods)
	// for which evaluation is computationally expensive and part of the computation can be shared
	// between the function and the derivative.
	//
	// The StoppingCondition is a function or function object
	// taking the evaluated point and function value as argument and returning true (stop iterating) or false (don't
	// stop iterating).
	//
	// The MatrixSolver is an object which can solve linear equations (i.e. decompose matrices).
	// This function is suitable for Eigen's solvers which support compute() and solve() methods.
	// A sensible option is usually Eigen::LDLT since the second derivative will be nonnegative-definite.
	//
	template< typename FunctionAndDerivativeEvaluator, typename StoppingCondition, typename MatrixSolver >
	typename FunctionAndDerivativeEvaluator::Vector find_root_by_newton_raphson(
		FunctionAndDerivativeEvaluator& evaluator,
		typename FunctionAndDerivativeEvaluator::Vector point,
		StoppingCondition& stopping_condition,
		MatrixSolver& solver
	)
	{
		typedef typename FunctionAndDerivativeEvaluator::Vector Vector ;
		typedef typename FunctionAndDerivativeEvaluator::Matrix Matrix ;
		
		evaluator.evaluate_at( point ) ;

		Matrix derivative_value ;

		Vector function_value = evaluator.get_value_of_function() ;
#if DEBUG_NEWTON_RAPHSON
		std::cerr << "NR: point = " << point << ", value = " << function_value << ".\n" ;
#endif
		while( !stopping_condition( point, function_value ) ) {
			// The Newton-Raphson rule comes from the observation that if
			// f( x + h ) = f( x ) + (D_x f) (h) + higher order terms
			// and if f( x + h ) = 0
			// then h must satisfy (D_x f) (h) = -f( x ) + higher order terms.
			// At each step we solve this and move to the point x + h.
			// If the function is linear, this will actually get us to the root.
			derivative_value = evaluator.get_value_of_first_derivative() ;
#if DEBUG_NEWTON_RAPHSON
			std::cerr << "NR: derivative_value = " << derivative_value << "\n" ;
#endif
			solver.compute( derivative_value ) ;
			point += solver.solve( -function_value ) ;
			evaluator.evaluate_at( point ) ;
#if DEBUG_NEWTON_RAPHSON
			std::cerr << "NR: point = " << point << ", value = " << function_value << ".\n" ;
#endif
			function_value = evaluator.get_value_of_function() ;
		}
		return point ;
	}
}

#endif
