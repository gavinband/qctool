
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef INTEGRATION_GANDERGAUTSCHIADAPTIVEINTEGRATOR_HPP
#define INTEGRATION_GANDERGAUTSCHIADAPTIVEINTEGRATOR_HPP

#include <vector>
#include <limits>
#include <cassert>
#include <iostream>
#include "Integrator.hpp"
#include "GanderGautschiAdaptiveIntegratorBase.hpp"

namespace integration {
	// class GanderGautschiAdaptiveIntegrator
	// This implements the Gauss-Kronod adaptive integration scheme from Gander-Gauschi,
	// "Adaptive quadrature - revisited" (1998)
	// Occasionally used notation below is as in the paper.
	class GanderGautschiAdaptiveIntegrator: public GanderGautschiAdaptiveIntegratorBase
	{
	public:
		GanderGautschiAdaptiveIntegrator( double desired_error ) ;

	private:
		double get_global_magnitude_estimate( double const xmin, double const xmax, double const global_integral_estimate ) const ;

		bool check_if_accuracy_met(
			std::vector< double > const& integral_estimates,
			double const global_magnitude_estimate
		) const ;

		std::vector< double > get_integral_estimates(
			std::vector< double > const& evaluations,
			double const h,
			std::size_t const number_of_estimates
		) const ;
		
		double combine_subintegrals( std::vector< double > const& subintegrals ) const ;
	} ;
}

#endif
