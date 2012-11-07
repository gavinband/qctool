
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "test_case.hpp"
#include "metro/FishersExactTest.hpp"

namespace  {
	void swap_rows( Eigen::Matrix2d& matrix ) {
		Eigen::Matrix2d P( 2, 2 ) ;
		P << 0, 1, 1, 0 ;
		matrix = P * matrix ;
	}
}

AUTO_TEST_CASE( test_fishers_exact_test_onesided ) {
	metro::FishersExactTest::Alternative const greater = metro::FishersExactTest::eGreater ;
	metro::FishersExactTest::Alternative const less = metro::FishersExactTest::eLess ;

	double const tolerance = 0.000001 ;
	Eigen::Matrix2d M ;
	M <<	1, 0,
			0, 1 ;
	BOOST_CHECK_EQUAL( metro::FishersExactTest( M ).get_pvalue( greater ), 0.5 ) ;
	BOOST_CHECK_EQUAL( metro::FishersExactTest( M ).get_pvalue( less ), 1.0 ) ;

	swap_rows( M ) ;
	BOOST_CHECK_EQUAL( metro::FishersExactTest( M ).get_pvalue( greater ), 1.0 ) ;
	BOOST_CHECK_EQUAL( metro::FishersExactTest( M ).get_pvalue( less ), 0.5 ) ;

	M <<	1, 1,
			1, 1 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( greater ), 0.833333333333333, tolerance ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( less ), 0.833333333333333, tolerance ) ;

	M <<	2, 0,
			1, 1 ;
	BOOST_CHECK_EQUAL( metro::FishersExactTest( M ).get_pvalue( greater ), 0.5 ) ;
	BOOST_CHECK_EQUAL( metro::FishersExactTest( M ).get_pvalue( less ), 1.0 ) ;

	swap_rows( M ) ;
	BOOST_CHECK_EQUAL( metro::FishersExactTest( M ).get_pvalue( greater ), 1.0 ) ;
	BOOST_CHECK_EQUAL( metro::FishersExactTest( M ).get_pvalue( less ), 0.5 ) ;

	M << 	15, 30,
			18, 27 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( greater ), 0.809109709011447, tolerance ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( less ), 0.3310533979, tolerance ) ;
	swap_rows( M ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( greater ), 0.3310533979, tolerance ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( less ), 0.809109709011447, tolerance ) ;

	M << 	15, 10,
			0, 10 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( greater ), 0.00100640923777742, tolerance ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( less ), 1, tolerance ) ;
	swap_rows( M ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( greater ), 1, tolerance ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( less ), 0.00100640923777742, tolerance ) ;

	M << 	15, 10,
			0, 20 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( greater ), 9.4783089312209e-06, tolerance ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( less ), 1, tolerance ) ;
	swap_rows( M ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( greater ), 1, tolerance ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( less ), 9.4783089312209e-06, tolerance ) ;

	M << 	15, 10,
			0, 30 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( greater ), 2.74692627172901e-07, tolerance ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( less ), 1, tolerance ) ;
	swap_rows( M ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( greater ), 1, tolerance ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( less ), 2.74692627172901e-07, tolerance ) ;
	
	M << 	1500, 1000,
			0, 1000 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( greater ), 5.614289492e-308, tolerance ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( less ), 1, tolerance ) ;
	swap_rows( M ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( greater ), 1, tolerance ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( less ), 5.614289492e-308, tolerance ) ;

	M <<	28, 2132,
			54, 1346 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( greater ), 1, tolerance ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( less ), 8.16627397356615e-07, tolerance ) ;
	swap_rows( M ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( greater ), 8.16627397356615e-07, tolerance ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( less ), 1, tolerance ) ;
}

AUTO_TEST_CASE( test_fishers_exact_test_twosided ) {
	metro::FishersExactTest::Alternative const twosided = metro::FishersExactTest::eTwoSided ;

	double const tolerance = 0.000001 ;
	Eigen::Matrix2d M ;
	M <<	1, 0,
			0, 1 ;
	BOOST_CHECK_EQUAL( metro::FishersExactTest( M ).get_pvalue( twosided ), 1.0 ) ;
	swap_rows(M) ;
	BOOST_CHECK_EQUAL( metro::FishersExactTest( M ).get_pvalue( twosided ), 1.0 ) ;

	M <<	1, 1,
				1, 1 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( twosided ), 1, tolerance ) ;
	swap_rows(M) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( twosided ), 1, tolerance ) ;

	M <<	2, 0,
			1, 1 ;
	BOOST_CHECK_EQUAL( metro::FishersExactTest( M ).get_pvalue( twosided ), 1 ) ;
	swap_rows(M) ;
	BOOST_CHECK_EQUAL( metro::FishersExactTest( M ).get_pvalue( twosided ), 1 ) ;

	M << 	15, 30,
			18, 27 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( twosided ), 0.6621067958, tolerance ) ;
	swap_rows(M) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( twosided ), 0.6621067958, tolerance ) ;

	M << 	15, 10,
			0, 10 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( twosided ), 0.001568035446, tolerance ) ;
	swap_rows(M) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( twosided ), 0.001568035446, tolerance ) ;

	M << 	15, 10,
			0, 20 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( twosided ), 1.233304071e-05, tolerance ) ;
	swap_rows(M) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( twosided ), 1.233304071e-05, tolerance ) ;

	M << 	15, 10,
			0, 30 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( twosided ), 2.746926272e-07, tolerance ) ;
	swap_rows(M) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( twosided ), 2.746926272e-07, tolerance ) ;

	M << 	1500, 1000,
			0, 1000 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( twosided ), 5.8781281e-308, tolerance ) ;
	swap_rows(M) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( twosided ), 5.8781281e-308, tolerance ) ;

	M <<	28, 2132,
			54, 1346 ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( twosided ), 1.00564929884086e-06, tolerance ) ;
	swap_rows( M ) ;
	BOOST_CHECK_CLOSE( metro::FishersExactTest( M ).get_pvalue( twosided ), 1.00564929884086e-06, tolerance ) ;
}

