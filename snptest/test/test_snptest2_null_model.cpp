
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include "Eigen/Eigen"
#include "test_case.hpp"
#include "snptest/case_control/NullModelLogLikelihood.hpp"
#include "snptest/FinitelySupportedFunctionSet.hpp"

AUTO_TEST_CASE( test_null_model_small_datasets )
{
	std::cerr << "Testing SNPTEST2 null model (small datasets)..." ;
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::MatrixXd Matrix ;
	using snptest::FinitelySupportedFunctionSet ;
	Vector levels( 3 ) ;
	levels << 0, 1, 2 ;

	{
		Vector const phenotypes = Vector::Zero( 1 ) ;
		Matrix genotypes = Matrix::Zero(1,3) ;
		genotypes(0,0) = 1.0 ; 
		snptest::case_control::NullModelLogLikelihood ll ;
		ll.set_phenotypes( phenotypes ).set_predictor_probs( genotypes, levels ) ;
		Vector parameters = Vector::Zero( 1 ) ;
		ll.evaluate_at( parameters ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( 0.5 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 1, -0.5 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), Matrix::Constant( 1, 1, -0.25) ) ;

		parameters( 0 ) = 1.0 ;
		double p0 = 1.0 / ( 1.0 + std::exp( 1.0 )) ;
		double p1 = std::exp(1.0) / ( 1.0 + std::exp( 1.0 )) ;
		ll.evaluate_at( parameters ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( p0) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 1, -p1) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), Matrix::Constant( 1, 1, -p1 * p0) ) ;
	}
	

	{
		Vector phenotypes( 2 ) ;
		phenotypes << 0.0, 1.0 ;
		Matrix genotypes = Matrix::Zero(2,3) ;
		genotypes(0,0) = 1.0 ; 
		genotypes(1,0) = 1.0 ; 

		snptest::case_control::NullModelLogLikelihood ll ;
		ll.set_phenotypes( phenotypes ).set_predictor_probs( genotypes, levels ) ;
		Vector parameters = Vector::Zero( 1 ) ;
		ll.evaluate_at( parameters ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( 0.5 ) + std::log( 0.5 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 1, 0.0 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), Matrix::Constant( 1, 1, -0.5) ) ;
		
		parameters( 0 ) = 1.0 ;
		double p0 = 1.0 / ( 1.0 + std::exp( 1.0 )) ;
		double p1 = std::exp(1.0) / ( 1.0 + std::exp( 1.0 )) ;
		ll.evaluate_at( parameters ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( p0 ) + std::log( p1) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 1, 1.0 - 2.0 * p1) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), Matrix::Constant( 1, 1, -2.0 * p1 * ( 1.0 - p1 )) ) ;
	}

	{
		Vector phenotypes( 2 ) ;
		phenotypes << 0.0, 1.0 ;
		Matrix genotypes = Matrix::Zero(2,3) ;
		genotypes(0,0) = 0.5 ; 
		genotypes(0,1) = 0.5 ; 
		genotypes(1,0) = 0.6 ; 
		genotypes(1,1) = 0.4 ; 

		snptest::case_control::NullModelLogLikelihood ll ;
		ll.set_phenotypes( phenotypes ).set_predictor_probs( genotypes, levels ) ;
		Vector parameters = Vector::Zero( 1 ) ;
		ll.evaluate_at( parameters ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( 0.5 ) + std::log( 0.5 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 1, 0.0 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), Matrix::Constant( 1, 1, -0.5) ) ;
		
		parameters( 0 ) = 1.0 ;
		double p0 = 1.0 / ( 1.0 + std::exp( 1.0 )) ;
		double p1 = std::exp(1.0) / ( 1.0 + std::exp( 1.0 )) ;
		ll.evaluate_at( parameters ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( p0 ) + std::log( p1) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 1, 1.0 - 2.0 * p1) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), Matrix::Constant( 1, 1, -2.0 * p1 * ( 1.0 - p1 )) ) ;
	}

	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_null_model_exclusions )
{
	std::cerr << "Testing SNPTEST2 null model (exclusions)..." ;
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::MatrixXd Matrix ;
	using snptest::FinitelySupportedFunctionSet ;
	Vector levels( 3 ) ;
	levels << 0, 1, 2 ;

    // Should get the same results by adding two indnividuals and excluding one.
	{
		Vector phenotypes( 2 ) ;
		phenotypes << 1.0, 0.0 ;
		Matrix genotypes = Matrix::Zero(2,3) ;
		genotypes(0,0) = 1.0 ; 
		genotypes(1,0) = 1.0 ; 
		std::vector< int > exclusions( 1, 0 ) ;

		snptest::case_control::NullModelLogLikelihood ll ;
		ll.set_phenotypes( phenotypes ).set_predictor_probs( genotypes, levels ).add_exclusions( exclusions ) ;
		ll.add_exclusions( exclusions ) ;
		Vector parameters = Vector::Zero( 1 ) ;
		ll.evaluate_at( parameters ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( 0.5 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 1, -0.5 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), Matrix::Constant( 1, 1, -0.25) ) ;

		parameters( 0 ) = 1.0 ;
		double p0 = 1.0 / ( 1.0 + std::exp( 1.0 )) ;
		double p1 = std::exp(1.0) / ( 1.0 + std::exp( 1.0 )) ;
		ll.evaluate_at( parameters ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( p0) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 1, -p1) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), Matrix::Constant( 1, 1, -p1 * p0) ) ;
	}

	// Should work the other way round too
	{
		Vector phenotypes( 2 ) ;
		phenotypes << 0.0, -1 ;
		Matrix genotypes = Matrix::Zero(2,3) ;
		genotypes(0,0) = 1.0 ; 
		genotypes(1,0) = 1.0 ; 
		std::vector< int > exclusions( 1, 1 ) ;

		snptest::case_control::NullModelLogLikelihood ll ;
		ll.set_phenotypes( phenotypes ).set_predictor_probs( genotypes, levels ).add_exclusions( exclusions ) ;
		
		Vector parameters = Vector::Zero( 1 ) ;
		ll.evaluate_at( parameters ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( 0.5 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 1, -0.5 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), Matrix::Constant( 1, 1, -0.25) ) ;

		parameters( 0 ) = 1.0 ;
		double p0 = 1.0 / ( 1.0 + std::exp( 1.0 )) ;
		double p1 = std::exp(1.0) / ( 1.0 + std::exp( 1.0 )) ;
		ll.evaluate_at( parameters ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( p0) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 1, -p1) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), Matrix::Constant( 1, 1, -p1 * p0) ) ;
	}

	// Try three smaples and two exclusions.
	{
		Vector phenotypes( 3 ) ;
		phenotypes << -1, 0.0, -1 ;
		Matrix genotypes = Matrix::Zero(3,3) ;
		genotypes(0,0) = 1.0 ; 
		genotypes(1,0) = 1.0 ; 
		genotypes(2,0) = 1.0 ; 
		std::vector< int > exclusions( 2 ) ;
        exclusions[0] = 0 ;
        exclusions[0] = 2 ;

		snptest::case_control::NullModelLogLikelihood ll ;
		ll.set_phenotypes( phenotypes ).set_predictor_probs( genotypes, levels ).add_exclusions( exclusions ) ;
		Vector parameters = Vector::Zero( 1 ) ;
		ll.evaluate_at( parameters ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( 0.5 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 1, -0.5 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), Matrix::Constant( 1, 1, -0.25) ) ;

		parameters( 0 ) = 1.0 ;
		double p0 = 1.0 / ( 1.0 + std::exp( 1.0 )) ;
		double p1 = std::exp(1.0) / ( 1.0 + std::exp( 1.0 )) ;
		ll.evaluate_at( parameters ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( p0) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 1, -p1) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), Matrix::Constant( 1, 1, -p1 * p0) ) ;
	}
}
