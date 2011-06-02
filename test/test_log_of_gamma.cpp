
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "gamma.hpp"
#include "floating_point_utils.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "test_case.hpp"

namespace {
	double naive_log_of_factorial( double x ) {
		assert( x > 0.0 ) ;
		if( x > 2.0 ) {
			return std::log( x ) + naive_log_of_factorial( x - 1.0 ) ;
		}
		else {
			return std::log( x ) ;
		}
	}

	void test_log_of_gamma( double tolerance ) {
		for( long x = 2; x < 10000; ++x ) {
			double a = naive_log_of_factorial( static_cast< double >( x - 1.0 )) ;
			double b = log_of_gamma_function( static_cast< double >( x )) ;
			assert( floats_are_equal( a, b, tolerance )) ;
		}
	}
}

AUTO_TEST_CASE( test_log_gamma ) {
	double tolerance = 0.0000000001 ;
	test_log_of_gamma( tolerance ) ;
}
