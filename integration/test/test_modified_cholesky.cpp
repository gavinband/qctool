
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
	Eigen::LDLT< Matrix > ldlt ;

	Matrix const identity = Matrix::Identity( 2, 2 ) ;

	Matrix test( 2, 2 ) ;
	test << 1, 2, 2, 1 ;

	{
		double offs[9] = { 0, 0.5, 0.9, 0.99, 1.0, 1.01, 1.05, 1.1, 2 } ;
		for( std::size_t i = 0; i < 9; ++i ) {
			test(0,1) = test(1,0) = offs[i] ;
			comp.compute( test ) ;
			BOOST_CHECK_EQUAL( comp.matrixP() * identity, identity ) ;
			// Diagonal should be strictly positive.
			for( int j = 0; j < 2; ++j ) {
				TEST_ASSERT( comp.vectorD()(j) >= std::numeric_limits< double >::epsilon() ) ;
			}

			std::cerr << "---------\n" ;
			std::cerr << "test_modified_cholesky_nonpositive_definite(): off-diagonal = " << test(0,1) << ",\n" ;
			std::cerr << "test_modified_cholesky_nonpositive_definite(): original matrix is:\n"
				<< test << ".\n" ;
			Matrix reconstructed = comp.matrixL() * Matrix( Vector( comp.vectorD() ).asDiagonal() ) * comp.matrixL().transpose() ;
			std::cerr << "test_modified_cholesky_nonpositive_definite(): reconstructed matrix is:\n"
				<< reconstructed << ".\n" ;
			
			// Modified cholesky is cholesky of a modified matrix with additions of +ve elements
			// to the diagonal.  Check changes are +ve
			BOOST_CHECK( (reconstructed.diagonal() - test.diagonal()).minCoeff() >= 0 ) ;
			// Check that diagonal changes are within error bound
			// (as in Fang & Leary (2006), formula (3)):
			double const beta = std::sqrt( std::max( 1.0, offs[i] / std::sqrt( test.rows() * test.rows() - 1 )) ) ;
			double const errorBound = (
				std::pow( offs[i] / beta + ( test.rows() - 1.0 ) * beta, 2 )
					+ 2 * ( test.diagonal().maxCoeff() + ( test.rows() - 1.0 ) * beta*beta )
					+ std::numeric_limits< double >::epsilon()
			) ;
			BOOST_CHECK( (reconstructed.diagonal() - test.diagonal()).maxCoeff() <= errorBound ) ;
			std::cerr << "test_modified_cholesky_nonpositive_definite(): bound is " << errorBound << ".\n" ;
			
			// Check that error is confined to diagonal
			{
				Matrix a = reconstructed.triangularView< Eigen::UnitLower >() ;
				Matrix b = test.triangularView< Eigen::UnitLower >() ;
				Matrix c = reconstructed.triangularView< Eigen::UnitUpper >() ;
				Matrix d = test.triangularView< Eigen::UnitUpper >() ;
				BOOST_CHECK( a == b ) ;
				BOOST_CHECK( c == d ) ;
			}
				
			if( test(0,1) < 1.0 ) {
				// Matrix is +ve definite, we should be able to reconstruct original matrix exactly
				BOOST_CHECK( (reconstructed.diagonal() - test.diagonal()).minCoeff() == 0 ) ;
				Matrix solved = comp.solve( identity ) ;
				BOOST_CHECK_EQUAL( solved * test, identity ) ;
			}
		}
	}
}


AUTO_TEST_CASE( test_modified_cholesky_special_test_cases ) {
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::MatrixXd Matrix ;
	Eigen::LDLT< Matrix > ldlt ;

	Matrix const identity = Matrix::Identity( 2, 2 ) ;

	Matrix test( 4, 4 ) ;
	
	// This is a problematic 4x4 matrix identified by Schnabel and Eskow.
	// It is discussed in section 6.2 of Fang & Leary (2006), who list
	// its smallest eigenvalue as -0.378 and the maximum error / min eigenvalue
	// as 2.733 for the GMW method.
	test <<
		 1890.3, -1705.6, -315.8,  3000.3,
		-1705.6,  1538.3,  284.9, -2706.6,
		 -315.8,   284.9,   52.5,  -501.2,
		 3000.3, -2706.6, -501.2,  4760.8 ;

	{
		integration::ModifiedCholesky< Matrix > comp ;
		comp.compute( test ) ;
		
		// Check diagonal is strictly positive.
		for( int j = 0; j < 2; ++j ) {
			TEST_ASSERT( comp.vectorD()(j) >= std::numeric_limits< double >::epsilon() ) ;
		}

		std::cerr << "---------\n" ;
		std::cerr << "test_modified_cholesky_special_test_cases(): off-diagonal = " << test(0,1) << ",\n" ;
		std::cerr << "test_modified_cholesky_special_test_cases(): original matrix is:\n"
			<< test << ".\n" ;
		Matrix reconstructed =
			comp.matrixL()
			* Matrix( Vector( comp.vectorD() ).asDiagonal() )
			* comp.matrixL().transpose() ;
		std::cerr << "test_modified_cholesky_special_test_cases(): reconstructed matrix (1) is:\n"
			<< reconstructed << ".\n" ;
		reconstructed = comp.matrixP().transpose() * reconstructed ;
		std::cerr << "test_modified_cholesky_special_test_cases(): reconstructed matrix (2) is:\n"
			<< reconstructed << ".\n" ;
		reconstructed = reconstructed * comp.matrixP() ;
		std::cerr << "test_modified_cholesky_special_test_cases(): reconstructed matrix (full) is:\n"
			<< reconstructed << ".\n" ;
		std::cerr << "test_modified_cholesky_special_test_cases(): permutation matrix is:\n"
			<< ( comp.matrixP() * Matrix::Identity( test.rows(), test.rows() ) ) << ".\n" ;
		
		// Modified cholesky is cholesky of a modified matrix with additions of +ve elements
		// to the diagonal.  Check changes are +ve
		BOOST_CHECK( (reconstructed.diagonal() - test.diagonal()).minCoeff() >= 0 ) ;
		// Check that diagonal changes are within error bound
		// (as in Fang & Leary (2006), formula (3)):
		double const beta = std::sqrt( std::max( 4760.8, 3000.3 / std::sqrt( test.rows() * test.rows() - 1 )) ) ;
		double const errorBound = (
			std::pow( 3000.3 / beta + ( test.rows() - 1.0 ) * beta, 2 )
				+ 2 * ( 4760.8 + ( test.rows() - 1.0 ) * beta*beta )
				+ std::numeric_limits< double >::epsilon()
		) ;
		double const fangLearyR2 = 2.7335 ; // Add 0.0005 for rounding
		Matrix const E = reconstructed - test ;
		std::cerr << "test_modified_cholesky_special_test_cases(): bound is " << errorBound << ", predicted r2 is: " << fangLearyR2 << "\n" ;
		std::cerr << "test_modified_cholesky_special_test_cases(): max error is " << E.maxCoeff() << ".\n" ;
		BOOST_CHECK( E.maxCoeff() <= errorBound ) ;
		BOOST_CHECK( E.maxCoeff() / 0.3780759 <= fangLearyR2 ) ;
		
		// Check that error is confined to diagonal
		{
			Matrix a = reconstructed.triangularView< Eigen::UnitLower >() ;
			Matrix b = test.triangularView< Eigen::UnitLower >() ;
			Matrix c = reconstructed.triangularView< Eigen::UnitUpper >() ;
			Matrix d = test.triangularView< Eigen::UnitUpper >() ;
			BOOST_CHECK( (a - b).array().abs().maxCoeff() < 1E-9 ) ;
			BOOST_CHECK( (c - d).array().abs().maxCoeff() < 1E-9 ) ;
		}
	}
}


