#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cassert>
#include <iostream>
#include "fputils/floating_point_utils.hpp"
#include "integration/Integrator.hpp"
#include "integration/GanderGautschiAdaptiveLogIntegrator.hpp"

namespace integration {
	GanderGautschiAdaptiveLogIntegrator::GanderGautschiAdaptiveLogIntegrator( double desired_error ):
		GanderGautschiAdaptiveIntegratorBase( desired_error )
	{}

	double GanderGautschiAdaptiveLogIntegrator::get_global_magnitude_estimate( double const xmin, double const xmax, double const global_integral_estimate ) const {
		if( global_integral_estimate > -std::numeric_limits< double >::infinity() ) {
			return global_integral_estimate ;
		}
		else {
			return std::log( xmax - xmin ) ;
		}
	}

	bool GanderGautschiAdaptiveLogIntegrator::check_if_accuracy_met(
		std::vector< double > const& integral_estimates,
		double const global_magnitude_estimate
	) const {
		// return( std::abs( integral_estimates[1] - integral_estimates[0] ) <= desired_error() * 0.1 * global_magnitude_estimate ) ;
			if( integral_estimates[0] < integral_estimates[1] ) {
				return ( fputils::log_diff_exp( integral_estimates[1], integral_estimates[0] ) <= std::log( desired_error() ) + global_magnitude_estimate ) ;
			}
			else {
				return ( fputils::log_diff_exp( integral_estimates[0], integral_estimates[1] ) <= std::log( desired_error() ) + global_magnitude_estimate ) ;
			}
	}

	std::vector< double > GanderGautschiAdaptiveLogIntegrator::get_integral_estimates(
		std::vector< double > const& evaluations,
		double const half_length_of_interval,
		std::size_t const number_of_estimates
	) const {
		assert( evaluations.size() % 2u == 1u ) ;
		std::vector< double > estimates( number_of_estimates ) ;
		for( std::size_t i = 0; i < number_of_estimates; ++i ) {
			std::vector< double > Is ;
			Is.reserve( evaluations.size() ) ;
			for( std::size_t j = 0; j < get_quadrature_rule( i ).size(); ++j ) {
				std::size_t const evaluation_point_index = get_quadrature_rule( i )[j].first ;
				double const weight = get_quadrature_rule( i )[ j ].second ;
				if( evaluation_point_index == ((evaluations.size() - 1) / 2 )) {
					Is.push_back( evaluations[ evaluation_point_index ] + std::log( weight ) ) ;
				}
				else {
					Is.push_back( evaluations[ evaluation_point_index ] + std::log( weight ) ) ;
					Is.push_back( evaluations[ evaluations.size() - 1 - evaluation_point_index ] + std::log( weight ) ) ;
				}
			}
			double I = fputils::log_sum_exp( Is.begin(), Is.end() ) ;
			// Above estimate is as though for integral from -1 to 1; adjust for actual length of interval.
			estimates[i] = I + std::log( half_length_of_interval ) ;
		}
		
		return estimates ;
	}
	
	double GanderGautschiAdaptiveLogIntegrator::combine_subintegrals( std::vector< double > const& subintegrals ) const {
		return fputils::log_sum_exp( subintegrals.begin(), subintegrals.end() ) ;
	}
}
