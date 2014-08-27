
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <cmath>
#include <utility>
#include <numeric>
#include "integration/Integrator.hpp"
#include "integration/GanderGautschiAdaptiveIntegrator.hpp"

namespace integration {
	GanderGautschiAdaptiveIntegrator::GanderGautschiAdaptiveIntegrator( double desired_error ):
		GanderGautschiAdaptiveIntegratorBase( desired_error )
	{}

	double GanderGautschiAdaptiveIntegrator::get_global_magnitude_estimate( double const xmin, double const xmax, double const global_integral_estimate ) const {
		double global_magnitude_estimate = std::abs( global_integral_estimate ) ;
		if( global_magnitude_estimate == 0.0 ) {
			global_magnitude_estimate = xmax - xmin ;
		}
		return global_magnitude_estimate ;
	}
	
	bool GanderGautschiAdaptiveIntegrator::check_if_accuracy_met(
		std::vector< double > const& integral_estimates,
		double const global_magnitude_estimate
	) const {
		// Gander and Gautschi compare the estimated error to the desired global error times the estimate of the
		// whole integrals' magnitude.  I add a factor of 0.1 because this more often leads to the actual accuracy
		// desired, (though it must also involve more function evaluations in general).
		return ( std::abs( integral_estimates[1] - integral_estimates[0] ) <= ( desired_error() * 0.1 * global_magnitude_estimate )) ;
	}

	std::vector< double > GanderGautschiAdaptiveIntegrator::get_integral_estimates(
		std::vector< double > const& evaluations,
		double const half_length_of_interval,
		std::size_t number_of_estimates
	) const {
		assert( evaluations.size() % 2u == 1u ) ;
		std::vector< double > estimates( number_of_estimates ) ;
		for( std::size_t i = 0; i < number_of_estimates; ++i ) {
			double I = 0.0 ;
			for( std::size_t j = 0; j < get_quadrature_rule( i ).size(); ++j ) {
				std::size_t const evaluation_point_index = get_quadrature_rule( i )[j].first ;
				double const weight = get_quadrature_rule( i )[ j ].second ;
				if( evaluation_point_index == ((evaluations.size() - 1) / 2 )) {
					I += evaluations[ evaluation_point_index ] * weight ;
				}
				else {
					I += ( evaluations[ evaluation_point_index ] + evaluations[ evaluations.size() - 1 - evaluation_point_index ] ) * weight ;
				}
			}
			// Above estimate is as though for integral from -1 to 1; adjust for actual length of interval.
			estimates[i] = I * half_length_of_interval ;
			if( I != I ) {
				throw IntegralIsNaNError() ;
			}
		}
		return estimates ;
	}
	
	double GanderGautschiAdaptiveIntegrator::combine_subintegrals( std::vector< double > const& subintegrals ) const {
		return std::accumulate( subintegrals.begin(), subintegrals.end(), 0.0 ) ;
	}
}

