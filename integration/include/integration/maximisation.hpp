
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_INTEGRATION_MAXIMISATION_HPP
#define QCTOOL_INTEGRATION_MAXIMISATION_HPP

#include "integration/NewtonRaphson.hpp"
#include "integration/Derivative.hpp"

namespace integration {
	template< typename SmoothFunction, typename Domain, typename StoppingCondition >
	Domain maximise_by_newton_raphson( SmoothFunction& f, Domain const& initial_point, StoppingCondition& stopping_condition ) {
		Derivative< SmoothFunction > Df = derivative( f ) ;
		return find_root_by_newton_raphson( Df, initial_point, stopping_condition ) ;
	}

	template< typename SmoothFunction, typename Domain >
	Domain maximise_by_newton_raphson( SmoothFunction& f, Domain const& initial_point, double tolerance = 0.0000000001, std::size_t max_iterations = 10000 ) {
		Derivative< SmoothFunction > Df = derivative( f ) ;
		return find_root_by_newton_raphson( Df, initial_point, tolerance, max_iterations ) ;
	}
}

#endif