#ifndef INTEGRATION_GANDERGAUTSCHIADAPTIVEINTEGRATORBASE_HPP
#define INTEGRATION_GANDERGAUTSCHIADAPTIVEINTEGRATORBASE_HPP

#include <vector>
#include <limits>
#include <cassert>
#include <iostream>
#include "Integrator.hpp"

namespace integration {
	// class GanderGautschiAdaptiveIntegratorBase
	// This implements the Gauss-Kronod adaptive integration scheme from Gander-Gauschi,
	// "Adaptive quadrature - revisited" (1998)
	// Occasionally used notation below is as in the paper.
	// Given the needed function evaluations, the actual evaluation of the integrals
	// is dispatched to a virtual method which must be implemented by a derived class.
	// The intent of this is that the integrand may be encoded somehow -- e.g. to work in log space.
	class GanderGautschiAdaptiveIntegratorBase: public Integrator
	{
	public:
		GanderGautschiAdaptiveIntegratorBase( double desired_error ) ;
		virtual ~GanderGautschiAdaptiveIntegratorBase() {}

		template< typename Function1D >
		double operator()( Function1D const& integrand, double xmin, double xmax ) const {
			assert( xmax > xmin ) ;
			return evaluate_recursively(
				integrand,
				xmin,
				xmax,
				integrand( xmin ),
				integrand( xmax )
			) ;
		}
	
	protected:

		typedef std::vector< std::pair< std::size_t, double > > QuadratureRule ;

		QuadratureRule const& get_quadrature_rule( std::size_t i ) const {
			assert( i < m_rules.size() ) ;
			return m_rules[i] ;
		}

		virtual std::vector< double > get_integral_estimates(
			std::vector< double > const& evaluations,
			double const h,
			std::size_t const number_of_estimates
		) const = 0 ;

		virtual double get_global_magnitude_estimate( double const xmin, double const xmax, double const global_integral_estimate ) const = 0 ;

		virtual bool check_if_accuracy_met(
			std::vector< double > const& integral_estimates,
			double const global_magnitude_estimate
		) const = 0 ;

		virtual double combine_subintegrals( std::vector< double > const& subintegrals ) const = 0 ;

	private:
		std::vector< double > m_evaluation_points ;
		std::vector< QuadratureRule > m_rules ;

		template< typename Function1D >
		double evaluate_recursively(
			Function1D const& integrand,
			double const xmin,
			double const xmax,
			double const evaluation_at_xmin,
			double const evaluation_at_xmax,
			double global_magnitude_estimate = std::numeric_limits< double >::quiet_NaN()
		) const {
			std::size_t const number_of_evaluations = ( m_evaluation_points.size() * 2 ) - 1 ;
			std::vector< double > points_at_which_to_evaluate( number_of_evaluations ) ;
			std::vector< double > evaluations( number_of_evaluations ) ;

			assert( xmax > xmin ) ;
			double const half_length_of_interval = 0.5 * ( xmax - xmin ) ;
			double const midpoint = 0.5 * ( xmax + xmin ) ;
			
			// Evaluate at the symmetrically paired points...
			points_at_which_to_evaluate.front() = xmin ;
			evaluations.front() = evaluation_at_xmin ;
			points_at_which_to_evaluate.back() = xmax ;
			evaluations.back() = evaluation_at_xmax ;

			// Only do the odd-numbered points if we are estimating the global integral...
			bool const perform_global_estimate = ( global_magnitude_estimate != global_magnitude_estimate ) ; // test for NaN.
			std::size_t increment = perform_global_estimate ? 1 : 2 ;
			for( std::size_t i = increment; i < ( m_evaluation_points.size() - 1 ); i += increment ) {
				points_at_which_to_evaluate[ i ] 								= midpoint - ( half_length_of_interval * m_evaluation_points[i] ) ;
				points_at_which_to_evaluate[ number_of_evaluations - 1 - i ] 	= midpoint + ( half_length_of_interval * m_evaluation_points[i] ) ;
				evaluations[ i ]												= integrand( points_at_which_to_evaluate[ i ]) ;
				evaluations[ number_of_evaluations - 1 - i ] 					= integrand( points_at_which_to_evaluate[ number_of_evaluations - 1 - i ] ) ;
			}
			
			// Evaluate at the centre point
			points_at_which_to_evaluate[ m_evaluation_points.size() - 1 ] = midpoint ;
			evaluations[ m_evaluation_points.size() - 1 ] = integrand( midpoint ) ;			
			
			// ...get the integral estimates...
			std::vector< double > integral_estimates = this->get_integral_estimates(
				evaluations,
				half_length_of_interval,
				perform_global_estimate ? m_rules.size() : m_rules.size() - 1
			) ;
			
			if( perform_global_estimate ) {
				global_magnitude_estimate = this->get_global_magnitude_estimate( xmin, xmax, integral_estimates.back() ) ;
			}
			
			// ...perform tolerance tests...
			if( this->check_if_accuracy_met( integral_estimates, global_magnitude_estimate ) ) {
				return integral_estimates[1] ;
			}
			else if( ( points_at_which_to_evaluate[2] <= xmin ) || ( points_at_which_to_evaluate[ number_of_evaluations - 3 ] >= xmax )) {
				throw IntervalTooSmallToSubdivideError( xmin, xmax ) ;
			}
			else {
				std::vector< double > subintegrals ;
				subintegrals.reserve( 6 ) ;
				for( std::size_t i = 0 ; i < points_at_which_to_evaluate.size() - 1 ; i += 2 ) {
					subintegrals.push_back(
						evaluate_recursively(
							integrand,
							points_at_which_to_evaluate[ i ],
							points_at_which_to_evaluate[ i + 2 ],
							evaluations[ i ],
							evaluations[ i + 2 ],
							global_magnitude_estimate
						)
					) ;
				}
				
				return this->combine_subintegrals( subintegrals ) ;
			}
		}
	} ;
}

#endif
