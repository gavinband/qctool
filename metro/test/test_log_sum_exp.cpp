
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <vector>
#include "test_case.hpp"
#include "metro/log_sum_exp.hpp"

namespace {
	template< typename F >
	bool floats_are_equal_to_within_epsilon( F a, F b, F epsilon ) {
		if( std::abs( a ) == std::numeric_limits< F >::infinity() || std::abs( b ) == std::numeric_limits< F >::infinity() )  {
			return a == b ;
		}
		else {
	    	return (std::abs( a - b ) < epsilon ) ;
		}
	}
}

AUTO_TEST_CASE( test_log_sum_exp ) {
	std::cerr << "test_log_sum_exp...\n" ;
	
	std::vector< std::vector< double > > test_data ;
	std::vector< double > expected_results ;
	
	double const infinity = std::numeric_limits< double >::infinity() ;
	
	// Empty sum -- result is log( 0 )
	test_data.push_back( std::vector< double >() ) ;
	expected_results.push_back( -infinity ) ;

	// Sum of one term -- result is that term.
	test_data.push_back( std::vector< double >( 1, 0.0 ) ) ;
	expected_results.push_back( 0.0 ) ;

	test_data.push_back( std::vector< double >( 1, -infinity ) ) ;
	expected_results.push_back( -infinity ) ;

	// Sum of two terms
	test_data.push_back( std::vector< double >( 2, 0.0 )) ;
	expected_results.push_back( std::log( 2.0 )) ;

	// Sum of -infinities
	test_data.push_back( std::vector< double >( 2, -infinity )) ;
	expected_results.push_back( -infinity ) ;

	test_data.push_back( std::vector< double >( 3, -infinity )) ;
	expected_results.push_back( -infinity ) ;

	// adding one or more -infinities to this should make no difference
	test_data.push_back( std::vector< double >( 3, 0.0 )) ;
	test_data.back()[0] =-infinity ;
	expected_results.push_back( std::log( 2.0 )) ;

	test_data.push_back( std::vector< double >( 3, 0.0 )) ;
	test_data.back()[1] =-infinity ;
	expected_results.push_back( std::log( 2.0 )) ;

	test_data.push_back( std::vector< double >( 3, 0.0 )) ;
	test_data.back()[2] =-infinity ;
	expected_results.push_back( std::log( 2.0 )) ;

	test_data.push_back( std::vector< double >( 4, 0.0 )) ;
	test_data.back()[0] =-infinity ;
	test_data.back()[1] =-infinity ;
	expected_results.push_back( std::log( 2.0 )) ;

	// With different values
	test_data.push_back( std::vector< double >( 2, 1.0 )) ;
	expected_results.push_back( std::log( std::exp( 1.0 ) + std::exp( 1.0 ))) ;

	assert( test_data.size() == expected_results.size() ) ;

	for( std::size_t i = 0; i < test_data.size(); ++i ) {
		/*
		std::cerr << "test data: " ;
		for( std::size_t j = 0; j < test_data[i].size(); ++j ) {
			std::cerr << test_data[i][j] << " " ;
		}
		std::cerr << "\n" ;
		*/
		double result = metro::log_sum_exp( test_data[i] ) ;
		TEST_ASSERT( floats_are_equal_to_within_epsilon( result, expected_results[i], 1e-12 )) ;
	}
	std::cerr << "success.\n" ;
}
