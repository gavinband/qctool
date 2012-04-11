
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_TEST_CASE_HPP
#define QCTOOL_TEST_CASE_HPP

#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>
#include "../../config.hpp"

#if HAVE_BOOST_UNIT_TEST_FRAMEWORK
	#include "boost/test/auto_unit_test.hpp"
	#include "boost/test/test_tools.hpp"
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
