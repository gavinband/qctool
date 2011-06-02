#define BOOST_TEST_MODULE snptest
#include "test_case.hpp"

bool floats_are_equal( double a, double b ) {
	double const epsilon = 0.0000000000001 ;
	
	bool result ;
	
	if( std::abs( a ) == std::numeric_limits< double >::infinity() || std::abs( b ) == std::numeric_limits< double >::infinity()) {
		result = ( a == b ) ;
	}
	else {
		result = std::abs( a - b ) < epsilon ;
	}
	if( !result ) {
		std::cerr << "floats not equal: " << a << " != " << b << ".\n" ;
	}
	return result ;
}
