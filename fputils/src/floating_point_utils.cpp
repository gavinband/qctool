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
