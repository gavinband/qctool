#ifndef GENFILE_TEST_CASE_HPP
#define GENFILE_TEST_CASE_HPP

#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>
#include "../../config.hpp"

#if HAVE_BOOST_UNIT_TEST
	#define BOOST_AUTO_TEST_MAIN
	#include "boost/test/auto_unit_test.hpp"
	#define AUTO_TEST_CASE( param ) BOOST_AUTO_TEST_CASE(param)
	#define TEST_ASSERT( param ) BOOST_ASSERT( param )
#else
	#define AUTO_TEST_CASE( param ) void param()
	#define TEST_ASSERT( param ) assert( param )
#endif	


#ifdef HAVE_BOOST_UNIT_TEST
#define AUTO_TEST_MAIN int dummy_func()
#else
#define AUTO_TEST_MAIN int main( int argc, char** argv )
#endif

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


#endif
