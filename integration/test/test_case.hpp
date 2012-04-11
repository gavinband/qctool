
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef STATFILE_TEST_CASE_HPP
#define STATFILE_TEST_CASE_HPP

#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>
#include "../../config.hpp"

#if HAVE_BOOST_UNIT_TEST_FRAMEWORK
	#define BOOST_TEST_MODULE integration
	#include "boost/test/unit_test.hpp"
	#define AUTO_TEST_CASE( param ) BOOST_AUTO_TEST_CASE(param)
	#define TEST_ASSERT( param ) BOOST_ASSERT( param )
	#define AUTO_TEST_MAIN void test_case_dummy_function_WILL_NOT_BE_CALLED()
#else
	#define AUTO_TEST_CASE( param ) void param()
	#define TEST_ASSERT( param ) assert( param )
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
