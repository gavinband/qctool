
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "fputils/log_space_matrix_multiply.hpp"
#include "test_case.hpp"

using fputils::Matrix ;

bool operator==( Matrix const& left, Matrix const& right ) {
	if( left.size1() != right.size1() || left.size2() != right.size2() ) {
		return false ;
	}
	
	for( std::size_t i = 0; i < left.size1(); ++i ) {
		for( std::size_t j = 0; j < left.size2(); ++j ) {
			if( left( i, j ) != right( i, j )) {
				return false ;
			}
		}
	}

	return true ;
}

bool operator!=( Matrix const& left, Matrix const& right ) {
	return !(left == right ) ;
}

AUTO_TEST_CASE( test_log_space_matrix_multiply ) {
	std::cerr << "test_log_space_matrix_multiply...\n" ;
	using fputils::ConstantMatrix ;
	typedef std::pair< Matrix, Matrix > P ;
	std::vector< P > multiplicands ;
	std::vector< Matrix > expected_results ;
	
	// multiply 1x1 matrices
	multiplicands.push_back(
		P(
			ConstantMatrix( 1, 1, 0.0 ),
			ConstantMatrix( 1, 1, 0.0 )
		)
	) ;
	
	expected_results.push_back( ConstantMatrix( 1, 1, 0.0 ) ) ;

	multiplicands.push_back(
		P(
			ConstantMatrix( 1, 1, 0.0 ),
			ConstantMatrix( 1, 1, 1.0 )
		)
	) ;
	
	expected_results.push_back( ConstantMatrix( 1, 1, 1.0 ) ) ;

	multiplicands.push_back(
		P(
			ConstantMatrix( 1, 1, 1.0 ),
			ConstantMatrix( 1, 1, 0.0 )
		)
	) ;
	
	expected_results.push_back( ConstantMatrix( 1, 1, 1.0 ) ) ;

	multiplicands.push_back(
		P(
			ConstantMatrix( 1, 1, 1.0 ),
			ConstantMatrix( 1, 1, 1.0 )
		)
	) ;
	
	expected_results.push_back( ConstantMatrix( 1, 1, 2.0 ) ) ;
	
	// multiply 1x2 by 2x1 matrices
	multiplicands.push_back(
		P(
			ConstantMatrix( 1, 2, 0.0 ),
			ConstantMatrix( 2, 1, 0.0 )
		)
	) ;
	
	expected_results.push_back( ConstantMatrix( 1, 1, std::log( 2.0 ))) ;

	multiplicands.push_back(
		P(
			ConstantMatrix( 1, 2, 0.0 ),
			ConstantMatrix( 2, 1, 1.0 )
		)
	) ;
	
	expected_results.push_back( ConstantMatrix( 1, 1, std::log( 2.0 ) + 1.0 )) ;

	multiplicands.push_back(
		P(
			ConstantMatrix( 1, 2, 1.0 ),
			ConstantMatrix( 2, 1, 0.0 )
		)
	) ;
	
	expected_results.push_back( ConstantMatrix( 1, 1, std::log( 2.0 ) + 1.0 )) ;

	multiplicands.push_back(
		P(
			ConstantMatrix( 1, 2, 1.0 ),
			ConstantMatrix( 2, 1, 1.0 )
		)
	) ;
	
	expected_results.push_back( ConstantMatrix( 1, 1, std::log( 2.0 ) + 2.0 )) ;
	
	// multiply 2x2 by 2x2 matrices
	multiplicands.push_back(
		P(
			ConstantMatrix( 2, 2, 0.0 ),
			ConstantMatrix( 2, 2, 0.0 )
		)
	) ;
	
	expected_results.push_back( ConstantMatrix( 2, 2, std::log( 2.0 ))) ;

	multiplicands.push_back(
		P(
			ConstantMatrix( 2, 2, 0.0 ),
			ConstantMatrix( 2, 2, 1.0 )
		)
	) ;
	
	expected_results.push_back( ConstantMatrix( 2, 2, std::log( 2.0 ) + 1.0 )) ;

	multiplicands.push_back(
		P(
			ConstantMatrix( 2, 2, 1.0 ),
			ConstantMatrix( 2, 2, 0.0 )
		)
	) ;
	
	expected_results.push_back( ConstantMatrix( 2, 2, std::log( 2.0 ) + 1.0 )) ;

	multiplicands.push_back(
		P(
			ConstantMatrix( 2, 2, 1.0 ),
			ConstantMatrix( 2, 2, 1.0 )
		)
	) ;
	
	expected_results.push_back( ConstantMatrix( 2, 2, std::log( 2.0 ) + 2.0 )) ;
	
	assert( multiplicands.size() == expected_results.size() ) ;
	for( std::size_t i = 0; i < multiplicands.size(); ++i ) {
		std::cerr << "Test case " << i+1 << " of " << multiplicands.size() << "...\n" ;
		Matrix m = fputils::log_space_matrix_multiply( multiplicands[i].first, multiplicands[i].second ) ;
		if( m != expected_results[i] ) {
			std::cerr << "Multiplying\n  " << multiplicands[i].first << "\n  by\n" << multiplicands[i].second << ",\n" ;
			std::cerr << "Expected:\n  " << expected_results[i] << "\n  but got\n" << m << ".\n" ;
		}
		TEST_ASSERT( m == expected_results[i] ) ;
	}
}
