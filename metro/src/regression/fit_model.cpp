
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <utility>
#include <Eigen/Core>
#include <limits>
#include <cassert>
#include "metro/regression/LogLikelihood.hpp"
#include "metro/ModifiedCholesky.hpp"
#include "metro/ModifiedNewtonRaphson.hpp"
#include "metro/regression/fit_model.hpp"
#include "metro/Snptest25StoppingCondition.hpp"

#define USE_MODIFIED_NR 1
#define DEBUG 1

namespace metro {
	namespace regression {
		Stepper::~Stepper() {}

		CholeskyStepper::CholeskyStepper( double tolerance, int max_iterations, Tracer tracer ):
			m_tolerance( tolerance ),
			m_max_iterations( max_iterations ),
			m_tracer( tracer ),
			m_iteration( -1 ),
			m_target_ll( -std::numeric_limits< double >::infinity() )
		{
			assert( tolerance > 0 ) ;
		}

		bool CholeskyStepper::diverged() const { return m_iteration >= m_max_iterations || m_target_ll != m_target_ll ; }
		std::size_t CholeskyStepper::number_of_iterations() const { return m_iteration ; }

		bool CholeskyStepper::step( Function& function, Vector const& point, Vector* result ) {
			++m_iteration ;
			function.evaluate_at( point, 2 ) ;
			double const ll = function.get_value_of_function() ;
			m_solver.compute( -function.get_value_of_second_derivative() ) ;
			double const derivativeOneNorm = function.get_value_of_first_derivative().array().abs().maxCoeff() ;
			Vector step = m_solver.solve( function.get_value_of_first_derivative() ) ;
			double directional_derivative = function.get_value_of_first_derivative().transpose() * step ;

			// We stop when the rescaled derivative is within tolerance
			// and we don't expect to improve the function value.
			bool converged = (
				(derivativeOneNorm < m_tolerance)
				&& ( ll > ( m_target_ll - m_tolerance ) )
				&& (directional_derivative < m_tolerance )
			) ;

			if( m_tracer ) {
				m_tracer( m_iteration, ll, m_target_ll, point, function.get_value_of_first_derivative(), step, converged ) ;
			}

			if(
				converged || (m_iteration) >= m_max_iterations 
			) {
				return false ;
			} else {
				(*result) = step ;
				m_target_ll = std::max( m_target_ll, ll ) ;
				return true ;
			}
		}

		std::pair< bool, int > fit_model(
			metro::SmoothFunction& ll,
			std::string const& model_name,
			Eigen::VectorXd const& starting_point,
			Stepper& stopping_condition,
			std::vector< std::string >* comments
		) {
			typedef metro::regression::Design::Matrix Matrix ;
			// Fit null model
			assert( starting_point.size() == ll.number_of_parameters() ) ;
			ll.evaluate_at( starting_point ) ;

			bool success = true ;
			int iterations = 0 ;
			// Compute negative definiteness condition number.
			Eigen::SelfAdjointEigenSolver< Matrix > eigenSolver( ll.get_value_of_second_derivative() ) ;
			if( eigenSolver.eigenvalues().maxCoeff() > 0 ) {
				success = false ;
				if( comments ) {
					comments->push_back( model_name + ":not-negative-definite" ) ;
				}
			} else {
	#if 1
			Eigen::VectorXd parameters = metro::find_maximum_by_modified_newton_raphson_with_line_search(
				ll, starting_point, stopping_condition
			) ;
	#else
				metro::Snptest25StoppingCondition< metro::regression::LogLikelihood > stopping_condition(
					ll,
					0.01,
					100,
					0
					//&get_ui_context().logger()
				) ;

				Eigen::ColPivHouseholderQR< Matrix > solver ;
				Eigen::VectorXd parameters = maximise_by_newton_raphson( ll, starting_point, stopping_condition, solver ) ;
	#endif

				if( stopping_condition.diverged() ) {
					if( comments ) {
						comments->push_back( model_name + ":model_fit_error:failed_to_converge_to_mle" ) ;
					}
				}

				success = !stopping_condition.diverged() ;
				iterations = stopping_condition.number_of_iterations() ;
			}

			return std::make_pair( success, iterations ) ;
		}
	}
}