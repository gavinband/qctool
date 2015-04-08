
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <boost/bind.hpp>
#include <boost/noncopyable.hpp>
#include "test_case.hpp"
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include "integration/ModifiedCholesky.hpp"

AUTO_TEST_CASE( test_modified_cholesky_positive_definite ) {
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::MatrixXd Matrix ;
    integration::ModifiedCholesky< Matrix > comp ;
	Eigen::LDLT< Matrix > ldlt ;
	Matrix const identity = Matrix::Identity( 2, 2 ) ;

	
	{
		Matrix const identity = Matrix::Identity( 1, 1 ) ;
		Matrix test( 1, 1 ) ;
		test << 1 ;
		comp.compute( test ) ;
		ldlt.compute( test ) ;
		BOOST_CHECK_EQUAL( comp.matrixP() * identity, ldlt.transpositionsP() * identity ) ;
		BOOST_CHECK_EQUAL( comp.vectorD(), ldlt.vectorD() ) ;
		BOOST_CHECK_EQUAL( comp.matrixL() * identity, ldlt.matrixL() * identity ) ;
		BOOST_CHECK_EQUAL( comp.solve( identity ) * test, identity ) ;
	}

	{
		Matrix test( 2, 2 ) ;
		test << 1, 0, 0, 1 ;
		comp.compute( test ) ;
		ldlt.compute( test ) ;
		BOOST_CHECK_EQUAL( comp.matrixP() * identity, ldlt.transpositionsP() * identity ) ;
		BOOST_CHECK_EQUAL( comp.vectorD(), ldlt.vectorD() ) ;
		BOOST_CHECK_EQUAL( comp.matrixL() * identity, ldlt.matrixL() * identity ) ;
		BOOST_CHECK_EQUAL( comp.solve( identity ) * test, identity ) ;
	}

	{
		// This is a non-positive definite matrix
		Matrix test( 2, 2 ) ;
		test << 1, 0.5, 0.5, 1 ;
		comp.compute( test ) ;
		ldlt.compute( test ) ;
		BOOST_CHECK_EQUAL( comp.matrixP() * identity, ldlt.transpositionsP() * identity ) ;
		BOOST_CHECK_EQUAL( comp.vectorD(), ldlt.vectorD() ) ;
		BOOST_CHECK_EQUAL( comp.matrixL() * identity, ldlt.matrixL() * identity ) ;
		BOOST_CHECK_EQUAL( comp.solve( identity ) * test, identity ) ;
	}
}

AUTO_TEST_CASE( test_modified_cholesky_nonpositive_definite ) {
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::MatrixXd Matrix ;
    integration::ModifiedCholesky< Matrix > comp ;

	Matrix const identity = Matrix::Identity( 2, 2 ) ;

	{
		// This is a non-positive definite matrix
		Matrix test( 2, 2 ) ;
		test << 1, 2, 2, 1 ;
		comp.compute( test ) ;
		BOOST_CHECK_EQUAL( comp.matrixP() * identity, identity ) ;
		// Diagonal should be at least as large as that of LDLT.
		for( int j = 0; j < 2; ++j ) {
			TEST_ASSERT( comp.vectorD()(j) >= std::numeric_limits< double >::epsilon() ) ;
		}

		std::cerr << "test_modified_cholesky_nonpositive_definite(): reconstructed matrix is:\n"
			<< comp.solve( identity ).inverse() << ".\n" ;
	}	
}


