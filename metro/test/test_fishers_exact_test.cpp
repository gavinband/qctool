
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "test_case.hpp"
#include "metro/FishersExactTest.hpp"

AUTO_TEST_CASE( test_fishers_exact_test ) {
	double const tolerance = 0.000001 ;
	Eigen::Matrix2d matrix ;
	matrix <<	1, 0,
				0, 1 ;
	BOOST_CHECK_EQUAL( metro::FishersExactTest( matrix ).get_pvalue(), 0.5 ) ;
	matrix <<	1, 1,
				1, 1 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( matrix ).get_pvalue(), 0.833333333333333, tolerance ) ;
	matrix <<	2, 0,
				1, 1 ;
	BOOST_CHECK_EQUAL( metro::FishersExactTest( matrix ).get_pvalue(), 0.5 ) ;
	matrix << 	15, 30,
				18, 27 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( matrix ).get_pvalue(), 0.809109709011447, tolerance ) ;
	matrix << 	15, 10,
				0, 10 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( matrix ).get_pvalue(), 0.00100640923777742, tolerance ) ;
	matrix << 	15, 10,
				0, 20 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( matrix ).get_pvalue(), 9.4783089312209e-06, tolerance ) ;
	matrix << 	15, 10,
				0, 30 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( matrix ).get_pvalue(), 2.74692627172901e-07, tolerance ) ;
	matrix << 	1500, 1000,
				0, 1000 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( matrix ).get_pvalue(), 5.614289492e-308, tolerance ) ;
}
