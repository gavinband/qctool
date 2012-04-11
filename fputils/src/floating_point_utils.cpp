
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "fputils/floating_point_utils.hpp"

namespace fputils {
	double log_sum_exp( double a, double b ) {
		if( b > a ) {
			std::swap( a, b ) ;
		}
		if( a == -std::numeric_limits< double >::infinity() ) {
			return a;
		}
		else {
			return a + std::log( 1.0 + std::exp( b - a )) ;
		}
	}
}
