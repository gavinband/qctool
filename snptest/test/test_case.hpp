#ifndef STATFILE_TEST_CASE_HPP
#define STATFILE_TEST_CASE_HPP

#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>
#include "../../config.hpp"

#if HAVE_BOOST_UNIT_TEST_FRAMEWORK
	#include "boost/test/unit_test.hpp"
	#define AUTO_TEST_CASE( param ) BOOST_AUTO_TEST_CASE(param)
	#define TEST_ASSERT( param ) BOOST_ASSERT( param )
#else
	#define AUTO_TEST_CASE( param ) void param()
	#define TEST_ASSERT( param ) assert( param )
	#define AUTO_TEST_MAIN int main( int argc, char** argv )
#endif	

bool floats_are_equal( double a, double b ) ;

#endif
