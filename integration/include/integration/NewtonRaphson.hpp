#ifndef INTEGRATION_NEWTON_RAPHSON_HPP
#define INTEGRATION_NEWTON_RAPHSON_HPP

#include <iomanip>
#include <boost/numeric/ublas/vector.hpp>
#include <Eigen/Dense>

namespace integration {
    namespace impl {
        template< typename Vector >
        struct VectorTraits
        {
        } ;

        template<>
        struct VectorTraits< Eigen::VectorXd >
        {
            typedef Eigen::MatrixXd Matrix ;
            static double squared_norm( Eigen::VectorXd const& v ) { return v.squaredNorm() ; }
        } ;
    }
	
	template< typename Vector, typename Function, typename Derivative >
	Vector find_root_by_newton_raphson(
		Function const& function,
		Derivative const& derivative,
		Vector point,
		double tolerance = 0.0000000001
	) {
        typedef typename impl::VectorTraits< Vector >::Matrix Matrix ;

	    // We compare against squared norm so square the tolerance here.
        //tolerance *= tolerance ;
		assert( tolerance > 0.0 ) ;
		
		Vector function_value = function( point ) ;
		double max = std::max( std::abs( function_value.minCoeff() ), std::abs( function_value.maxCoeff() ) ) ;
		if( max >= tolerance ) {
			// The Newton-Raphson rule comes from the observation that if
			// f( x + h ) = f( x ) + (D_x f) (h) + higher order terms
			// and if f( x + h ) = 0
			// then h must satisfy (D_x f) (h) = -f( x ) + higher order terms.
			// At each step we solve this and move to the point x + h.
			// If the function is linear, this will actually get us to the root.
			Eigen::ColPivHouseholderQR< Matrix > decomposer ;
			do {
				decomposer.compute( derivative( point ) ) ;
                // The following line does not work with Eigen beta 1
				//point += decomposer.solve( -function_value ) ; // 
				point = point + decomposer.solve( -function_value ) ;
                function_value = function( point ) ;
				max = std::max( std::abs( function_value.minCoeff() ), std::abs( function_value.maxCoeff() ) ) ;				
				//std::cerr << "NR: point = " << point << ".\n" ;
				//std::cerr << "NR: tolerance = " << tolerance << ", value = " << function_value << ", max coeff = " << max << ".\n" ;
			}
            while( max > tolerance ) ;
		}
		return point ;
	}

	template< typename FunctionAndDerivativeEvaluator >
	typename FunctionAndDerivativeEvaluator::Point find_root_by_newton_raphson(
		FunctionAndDerivativeEvaluator& evaluator,
		typename FunctionAndDerivativeEvaluator::Point point,
		double tolerance = 0.0000000001
	)
	// The version of Newton-Raphson taking a seperate function and derivative argument
	// has the disadvantage that any calculations that are common between function
	// and derivative evaluations, cannot easily be shared.  This version allows
	// this work to be shared by using a single object to compute both function and derivative.
	// The evaluator must expose an Evaluation typedef.  This is an object with two methods,
	// get_value_of_function() and get_value_of_derivative().
	{
		typedef typename FunctionAndDerivativeEvaluator::Point Point ;
		typedef typename FunctionAndDerivativeEvaluator::Vector Vector ;
		typedef typename FunctionAndDerivativeEvaluator::Matrix Matrix ;
		
	    // We compare against squared norm so square the tolerance here.
        //tolerance *= tolerance ;
		assert( tolerance > 0.0 ) ;
		evaluator.evaluate_at( point ) ;
		Vector function_value = evaluator.get_value_of_function() ;
		Matrix derivative_value ;
		std::cerr << std::setprecision( 6 ) ;
		std::cerr << "find_root_by_newton_raphson(): point = " << point << ", function = " << function_value ;
		
		double max = std::max( std::abs( function_value.minCoeff() ), std::abs( function_value.maxCoeff() ) ) ;
		if( max >= tolerance ) {
			// The Newton-Raphson rule comes from the observation that if
			// f( x + h ) = f( x ) + (D_x f) (h) + higher order terms
			// and if f( x + h ) = 0
			// then h must satisfy (D_x f) (h) = -f( x ) + higher order terms.
			// At each step we solve this and move to the point x + h.
			// If the function is linear, this will actually get us to the root.
			Eigen::ColPivHouseholderQR< Matrix > solver ;
			do {
				derivative_value = evaluator.get_value_of_first_derivative() ;
				solver.compute( derivative_value ) ;
				std::cerr << ", derivative value = {" << derivative_value << "}.\n" ;
				std::cerr << ", addition = {" << solver.solve( -function_value ) << "}.\n" ;
                // The following line does not work with Eigen beta 1
				//point += solver.solve( -function_value ) ; // 
				point = point + solver.solve( -function_value ) ;
				evaluator.evaluate_at( point ) ;
				function_value = evaluator.get_value_of_function() ;
				max = std::max( std::abs( function_value.minCoeff() ), std::abs( function_value.maxCoeff() ) ) ;				
				//std::cerr << "NR: point = " << point << ".\n" ;
				//std::cerr << "NR: tolerance = " << tolerance << ", value = " << function_value << ", max coeff = " << max << ".\n" ;
				std::cerr << "find_root_by_newton_raphson(): point = " << point
					<< ", function = " << function_value ;
			}
            while( max > tolerance ) ;
			std::cerr << ".\n" ;
		}
		return point ;
	}

    // Specialise for 1d, where Vector == double
	template< typename Function, typename Derivative >
	double find_root_by_newton_raphson (
		Function const& function,
		Derivative const& derivative,
		double point,
		double const tolerance
	) {
		double function_value = function( point ) ;
		if( std::abs( function_value ) >= tolerance ) {
			// The Newton-Raphson rule comes from the observation that if
			// f( x + h ) = f( x ) + (D_x f) (h) + higher order terms
			// and if f( x + h ) = 0
			// then h must satisfy (D_x f) (h) = -f( x ) + higher order terms.
			// i.e. h ~ - f(x) / D_x f since we are in the scalar case.
			// At each step we solve this and move to the point x + h.
			// If the function is linear or quadratic, this will actually get us there.
			do {
				point -= function_value / derivative( point ) ;
                function_value = function( point ) ;
			}
            while( std::abs( function_value ) >= tolerance ) ;
		}
		return point ;
	}
}

#endif
