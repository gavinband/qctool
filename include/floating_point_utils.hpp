
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef __GTOOL_FLOATING_POINT_UTILS__
#define __GTOOL_FLOATING_POINT_UTILS__

#include <cmath>

template< typename F >
bool floats_are_equal( F a, F b, F tolerance ) {
	//See http://www.boost.org/doc/libs/1_34_1/libs/test/doc/components/test_tools/floating_point_comparison.html
    F difference = std::abs( a - b ) ;
    return (difference <= tolerance * std::abs( a ))
        && (difference <= tolerance * std::abs( b )) ;
}

template< typename F >
bool floats_are_equal_to_within_epsilon( F a, F b, F epsilon ) {
    return (std::abs( a - b ) < epsilon ) ;
}


template< typename F >
bool check_if_less_than_with_tolerance( F a, F b, F tolerance ) {
	return a < (b + tolerance) ;
}

template< typename F >
bool check_if_even( F a, F tolerance = 0.0 ) {
	return floats_are_equal( std::fmod( a, 2.0 ), 0.0, tolerance ) ;
}

template< typename F >
bool check_if_odd( F a, F tolerance = 0.0 ) {
	return floats_are_equal( std::fmod( a, 2.0 ), 1.0, tolerance ) ;
}

template< typename F >
F round_to_nearest_integer( F a ) {
	F result = std::floor( a ) ;
	if( a - result > 0.5 ) {
		result += 1.0 ;
	}
	return result ;
}

#endif

