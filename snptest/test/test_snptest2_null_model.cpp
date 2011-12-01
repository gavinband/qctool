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
		snptest::case_control::NullModelLogLikelihood ll( phenotypes, FinitelySupportedFunctionSet( levels, genotypes ) ) ;
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
	
	// Should get the same results by adding two indnividuals and excluding one.
	{
		Vector phenotypes( 2 ) ;
		phenotypes << 1.0, 0.0 ;
		Matrix genotypes = Matrix::Zero(2,3) ;
		genotypes(0,0) = 1.0 ; 
		genotypes(1,0) = 1.0 ; 
		std::vector< std::size_t > exclusions( 1, 0 ) ;

		snptest::case_control::NullModelLogLikelihood ll( phenotypes, FinitelySupportedFunctionSet( levels, genotypes ), Matrix(), exclusions ) ;
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

	// Should be true even if we have NaNs
	{
		Vector phenotypes( 2 ) ;
		phenotypes << 0.0, std::numeric_limits< double >::quiet_NaN() ;
		Matrix genotypes = Matrix::Zero(2,3) ;
		genotypes(0,0) = 1.0 ; 
		genotypes(1,0) = 1.0 ; 
		std::vector< std::size_t > exclusions( 1, 1 ) ;

		snptest::case_control::NullModelLogLikelihood ll( phenotypes, FinitelySupportedFunctionSet( levels, genotypes ), Matrix(), exclusions ) ;
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

		snptest::case_control::NullModelLogLikelihood ll( phenotypes, FinitelySupportedFunctionSet( levels, genotypes ) ) ;
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

