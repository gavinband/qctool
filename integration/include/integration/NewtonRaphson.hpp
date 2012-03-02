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
		double tolerance
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
				// std::cerr << "NR: point = " << point << ".\n" ;
				// std::cerr << "NR: tolerance = " << tolerance << ", value = " << function_value << ", max coeff = " << max << ".\n" ;
			}
            while( max > tolerance ) ;
		}
		return point ;
	}

	template< typename FunctionAndDerivativeEvaluator, typename StoppingCondition >
	typename FunctionAndDerivativeEvaluator::Vector find_root_by_newton_raphson(
		FunctionAndDerivativeEvaluator& evaluator,
		typename FunctionAndDerivativeEvaluator::Vector point,
		StoppingCondition& stopping_condition
	)
	// The version of Newton-Raphson taking a seperate function and derivative argument
	// has the disadvantage that any calculations that are common between function
	// and derivative evaluations, cannot easily be shared.  This version allows
	// this work to be shared by using a single object to compute both function and derivative.
	// The evaluator must expose an Evaluation typedef.  This is an object with two methods,
	// get_value_of_function() and get_value_of_derivative().
	{
		typedef typename FunctionAndDerivativeEvaluator::Vector Vector ;
		typedef typename FunctionAndDerivativeEvaluator::Matrix Matrix ;
		
		evaluator.evaluate_at( point ) ;

		Matrix derivative_value ;
		Eigen::ColPivHouseholderQR< Matrix > solver ;

		Vector function_value = evaluator.get_value_of_function() ;
		while( !stopping_condition( function_value ) ) {
			// The Newton-Raphson rule comes from the observation that if
			// f( x + h ) = f( x ) + (D_x f) (h) + higher order terms
			// and if f( x + h ) = 0
			// then h must satisfy (D_x f) (h) = -f( x ) + higher order terms.
			// At each step we solve this and move to the point x + h.
			// If the function is linear, this will actually get us to the root.
			derivative_value = evaluator.get_value_of_first_derivative() ;
			solver.compute( derivative_value ) ;
			point += solver.solve( -function_value ) ;
			evaluator.evaluate_at( point ) ;
				// std::cerr << "NR: point = " << point << ".\n" ;
				// std::cerr << "NR: tolerance = " << tolerance << ", value = " << function_value << ", max coeff = " << max << ".\n" ;
			function_value = evaluator.get_value_of_function() ;
		}
		return point ;
	}

	namespace impl {
		template< typename FunctionAndDerivativeEvaluator >
		struct FunctionNearZeroStoppingCondition
		{
			typedef typename FunctionAndDerivativeEvaluator::Vector Vector ;
			typedef typename FunctionAndDerivativeEvaluator::Matrix Matrix ;
			FunctionNearZeroStoppingCondition( double tolerance, std::size_t max_iterations ): m_tolerance( tolerance ), m_max_iterations( 10000 ), m_iteration( 0 ) {}
			bool operator()(
				Vector const& value_of_function
			) {
				std::cerr << "iteration " << m_iteration << ": value is " << std::resetiosflags( std::ios::floatfield ) << value_of_function << ".\n" ;
				return
					( ++m_iteration > m_max_iterations )
					||
					( std::max( std::abs( value_of_function.minCoeff() ), std::abs( value_of_function.maxCoeff() ) ) < m_tolerance )
				;
			}

		private:
			double const m_tolerance ;
			std::size_t const m_max_iterations ;
			std::size_t m_iteration ;
		} ;
	}

	template< typename FunctionAndDerivativeEvaluator >
	typename FunctionAndDerivativeEvaluator::Vector find_root_by_newton_raphson(
		FunctionAndDerivativeEvaluator& evaluator,
		typename FunctionAndDerivativeEvaluator::Vector point,
		double tolerance = 0.0000000001,
		std::size_t max_iterations = 10000
	) {
		assert( tolerance > 0.0 ) ;
		impl::FunctionNearZeroStoppingCondition< FunctionAndDerivativeEvaluator > stopping_condition( tolerance, max_iterations ) ;
		return find_root_by_newton_raphson(
			evaluator,
			point,
			stopping_condition
		) ;
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
