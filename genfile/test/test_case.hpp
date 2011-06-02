#ifndef GENFILE_TEST_CASE_HPP
#define GENFILE_TEST_CASE_HPP

#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>
#include "../../config.hpp"

#if HAVE_BOOST_UNIT_TEST_FRAMEWORK
	#include "boost/test/auto_unit_test.hpp"
	#define AUTO_TEST_CASE( param ) BOOST_AUTO_TEST_CASE(param)
	#define TEST_ASSERT( param ) BOOST_ASSERT( param )
	#define AUTO_TEST_MAIN namespace { void test_case_dummy_function_WILL_NOT_BE_CALLED() ; } void test_case_dummy_function_WILL_NOT_BE_CALLED() 
#else
	#define AUTO_TEST_CASE( param ) void param()
	#define TEST_ASSERT( param ) assert( param )
	#define AUTO_TEST_MAIN int main( int argc, char** argv )
#endif	

namespace data {
	std::string to_hex( std::string const& str ) ;
}

#endif
