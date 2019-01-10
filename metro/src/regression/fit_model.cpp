
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

			// Compute derivative normalised to standard Gaussian space.
			// Theory: if -I = LLᵗ then Σ=(Lᵗ)⁻¹L⁻¹.
			// We assume the log-likelihood looks like this:
			// f(x) = s(Lᵗ(x-x₀)) + O((x-x₀)³)
			// near the true maximum x₀, where s is the standard multivariate normal
			// log-density.  Then
			// f'(x) = -L Lᵗ (x-x₀) + O((x-x₀)²)   (expressed as a column vector).
			// We put a convergence condition on the corresponding derivative in s-space, given by
			// z = L⁻¹ f'(x) = -Lᵗ(x-x₀) + O((x-x₀)²)
			// This has the nice interpretation of being invariant to changes in scale of the
			// domain or range of the function, and also of being interpretable as a distance
			// to the maximum in the uncorrelated-variable space that is the domain of s.

			// Note further that the next step is (Lᵗ)⁻¹ (L⁻¹ f'(x)) so we have to
			// perform this computation anyway.  Conceptually a condition on the magnitude
			// of z is half-way between a condition on the magnitude of f' and on the magnitude
			// of the next jump.  A condition on z is independent of the scale of the function.

			// What could go wrong with this scheme?
			// If the function is not well-approximated by a quadratic then 
			// If the derivative is large but L has large entries it could be that 
			// z is small.

			//
			// Could this scheme go wrong away from the maximum?  Well it might stop early if
			// L⁻¹ has small values.
			// TODO: use special form of regression problem to show this is avoided.
			Vector step = function.get_value_of_first_derivative() ;
			double const derivativeOneNorm = step.array().abs().maxCoeff() ;
			m_solver.matrixL().solveInPlace( step ) ;
			//double const rescaledDerivativeOneNorm = step.array().abs().maxCoeff() ;

			// Compute the next step and the directional derivative
			m_solver.matrixL().transpose().solveInPlace( step ) ;
			double directional_derivative = function.get_value_of_first_derivative().transpose() * step ;

			// We stop when the rescaled derivative is within tolerance
			// and we don't expect to improve the function value.
			bool converged = (
				(derivativeOneNorm < m_tolerance)
				&& ( ll > ( m_target_ll - m_tolerance ) )
				&& (directional_derivative < m_tolerance )
			) ;

			if( m_tracer ) {
				m_tracer( m_iteration, m_target_ll, point, function.get_value_of_first_derivative(), step, converged ) ;
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
					options().get< double >( "-tolerance" ),
					options().get< std::size_t >( "-max-iterations" ),
					(appcontext::OstreamTee*)(0)
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