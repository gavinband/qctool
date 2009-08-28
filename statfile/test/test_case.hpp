#ifndef STATFILE_TEST_CASE_HPP
#define STATFILE_TEST_CASE_HPP

#include <cassert>
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
#define AUTO_TEST_MAIN void dummy_func()
#else
#define AUTO_TEST_MAIN int main( int argc, char** argv )
#endif

#endif
