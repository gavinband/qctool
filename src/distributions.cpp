/*
 * Note: classes in this file are designed to match boost math toolkit naming scheme.
 * There are two reasons for this choice:
 * 1. At some point, these may be replaced with the ones from boost.
 * 2. Even if not, we may as well follow an already-defined, well-thought-out api such as the one in boost.
 */

 
#include <cassert>
#include "distributions.hpp"
#include <cstddef>

namespace impl {
	double lookup_table[3][10] = {
		/* p-value  0.05    0.01    0.001 */
		/* 1 */		3.84, 	6.64, 	10.83,
		/* 2 */		5.99, 	9.21, 	13.82,
		/* 3 */		7.82, 	11.35, 	16.27,
		/* 4 */		9.49, 	13.28, 	18.47,
		/* 5 */		11.07, 	15.09, 	20.52,
		/* 6 */		12.59, 	16.81, 	22.46,
		/* 7 */		14.07, 	18.48, 	24.32,
		/* 8 */		15.51, 	20.09, 	26.13,
		/* 9 */		16.92, 	21.67, 	27.88,
		/* 10 */	18.31, 	23.21, 	29.59
	} ;

	double perform_lookup( double x, unsigned long degrees_of_freedom ) {
		std::size_t ypos = (degrees_of_freedom - 1) * 3;
		if( x >= lookup_table[ypos][2] ) {
			return 0.001 ;
		}
		else if( x >= lookup_table[ypos][1] ) {
			return 0.01 ;
		}
		else if( x >= lookup_table[ypos][0] ) {
			return 0.05 ;
		}
		else {
			return 1.0 ;
		}
	}
 }
 
// Return a value of x so that the probability of a chi-squared distributed variable
// lying to the right of x is at most the given p-value.
template<>
double quantile< complement_type< chi_squared_distribution > >( complement_type< chi_squared_distribution > const& distribution, double x ) {
	return impl::perform_lookup( x, distribution.degrees_of_freedom() ) ;
}
