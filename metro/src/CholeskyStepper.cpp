
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
#include "metro/CholeskyStepper.hpp"

#define DEBUG 1

namespace metro {
	CholeskyStepper::CholeskyStepper( double tolerance, int max_iterations, Tracer tracer ):
		m_tolerance( tolerance ),
		m_max_iterations( max_iterations ),
		m_tracer( tracer ),
		m_iteration( -1 ),
		m_target_ll( -std::numeric_limits< double >::infinity() )
	{
		assert( tolerance > 0 ) ;
	}
	
	bool CholeskyStepper::step( Function& function, Vector const& point, Vector* result ) {
		++m_iteration ;
		function.evaluate_at( point, 2 ) ;
		double const ll = function.get_value_of_function() ;
		m_solver.compute( -function.get_value_of_second_derivative() ) ;

		// Different possible stopping conditions are possible.
		// One is to stop when the predicted improvement in function
		// value is small.  This amount is encoded by the directional derivative.
		// However, for statistical model fitting we would like a rule that is invariant to
		// various types of rescaling the function.  For example
		// - scaling parameters (which scales the derivative)
		// - multiplying the function by a constant (which is a bit like adding more data)
		// A rule like directional_derivative < tolerance does not have these properties.

		// Instead, we compute the derivative normalised to 'standard Gaussian' space.
		// Theory: if the function 2nd derivative H = LLᵗ then Σ=(Lᵗ)⁻¹L⁻¹.
		// We assume the log-likelihood looks like this:
		// f(x) = s(Lᵗ(x-x₀)) + O((x-x₀)³)
		// near the true maximum x₀, where s is the standard multivariate normal
		// log-density.  Then
		// f'(x) = -L Lᵗ (x-x₀) + O((x-x₀)²)   (expressed as a column vector).
		// We put a convergence condition on the corresponding derivative in s-space, given by
		// z = L⁻¹ f'(x) = -Lᵗ(x-x₀) + O((x-x₀)²)
		// This has the nice interpretation of being interpretable as a distance
		// to the maximum in the uncorrelated-variable space that is the domain of s.

		double const derivativeOneNorm = function.get_value_of_first_derivative().array().abs().maxCoeff() ;
		Vector step = m_solver.solve( function.get_value_of_first_derivative() ) ;
		Vector const sqrtStep = m_solver.halfSolve( function.get_value_of_first_derivative() ) ;
//		double directional_derivative = function.get_value_of_first_derivative().transpose() * step ;

		bool converged = (
			(derivativeOneNorm < m_tolerance)
			&& ( ll >= ( m_target_ll - m_tolerance ) )
			&& (sqrtStep.array().abs().maxCoeff() < m_tolerance )
			//&& (directional_derivative < m_tolerance )
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
	
	bool CholeskyStepper::diverged() const {
		return m_iteration >= m_max_iterations || m_target_ll != m_target_ll ;
	}

	std::size_t CholeskyStepper::number_of_iterations() const {
		return m_iteration ;
	}

	void CholeskyStepper::reset() {
		m_iteration = -1 ;
		m_target_ll = -std::numeric_limits< double >::infinity() ;
	}
}
