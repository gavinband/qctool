#ifndef QCTOOL_INTEGRATION_MAXIMISATION_HPP
#define QCTOOL_INTEGRATION_MAXIMISATION_HPP

#include "integration/NewtonRaphson.hpp"
#include "integration/Derivative.hpp"

namespace integration {
	template< typename SmoothFunction, typename Domain >
	Domain maximise_by_newton_raphson( SmoothFunction& f, Domain const& initial_point, double tolerance ) {
		Derivative< SmoothFunction > Df = derivative( f ) ;
		return find_root_by_newton_raphson( Df, initial_point, tolerance ) ;
	}
}

#endif