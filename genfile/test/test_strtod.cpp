
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "test_case.hpp"
#include <string>
#include <cstdlib>
#include "genfile/string_utils/slice.hpp"
#include "genfile/string_utils/strtod.hpp"

namespace {
	void test_it( std::string const& value ) {
		double a = 0.0 ;
		double b = 0.0 ;
		bool error = false ;
		try {
			char* endptr ;
			char const* c_str = value.c_str() ;
			a = std::strtod( c_str, &endptr ) ;
			if( endptr != c_str + value.size() ) {
				error = true ;
			}
		}
		catch( std::exception const& e ) {
			error = true ;
		}

		try {
			b = genfile::string_utils::strtod( genfile::string_utils::slice( value ) ) ;
			TEST_ASSERT( !error ) ;
		}
		catch( std::exception const& e ) {
			TEST_ASSERT( error ) ;
		}
		
		if( (( a == a ) && ( b == b )) && ( a != b ) ) {
			std::cerr << "Comparison false.  a = " << a << ", b = " << b << ".\n" ;
		}
		
		TEST_ASSERT( (( a != a ) && ( b != b )) || ( a == b )) ;
	}
}

AUTO_TEST_CASE( test_strtod ) {
	test_it( "0" ) ;
	test_it( "nan" ) ;
	test_it( "naN" ) ;
	test_it( "nAn" ) ;
	test_it( "nAN" ) ;
	test_it( "Nan" ) ;
	test_it( "NaN" ) ;
	test_it( "NAn" ) ;
	test_it( "NAN" ) ;

	test_it( "inf" ) ;
	test_it( "inF" ) ;
	test_it( "iNf" ) ;
	test_it( "iNF" ) ;
	test_it( "Inf" ) ;
	test_it( "InF" ) ;
	test_it( "INf" ) ;
	test_it( "INF" ) ;

	test_it( "NA" ) ;

	for( double x = -100.0; x < 100.0; x += 0.0001 ) {
		std::ostringstream ostr ;
		ostr << x ;
		test_it( ostr.str() ) ;
	}
	
	// Here are some values known to cause problems for some implementations.
	test_it( "225073858507201e-308" ) ; // http://code.google.com/p/mochiweb/issues/detail?id=59#c0	
	test_it( "1.15507e-173" ) ; // found using test_strtod_random
}


AUTO_TEST_CASE( test_strtod_random ) {
	assert( sizeof( int ) == 4 ) ;
	assert( sizeof( double ) == 8 ) ;
	int D[2] ;
	for( std::size_t i = 0; i < 1000000; ++i ) {
		D[0] = rand() ;
		D[1] = rand() ;
		double x = *(reinterpret_cast< double* >( &D[0] )) ;
		std::ostringstream ostr ;
		ostr << x ;
		
		test_it( ostr.str() ) ;
	}
}
