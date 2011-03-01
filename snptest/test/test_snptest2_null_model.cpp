#include <iostream>
#include <iomanip>
#include "test_case.hpp"
#include "snptest/SNPTEST2NullModel.hpp"
#include "Eigen/Eigen"

AUTO_TEST_CASE( test_null_model_small_datasets )
{
	std::cerr << "Testing SNPTEST2 null model (small datasets)..." ;
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::MatrixXd Matrix ;
	
	{
		Vector const phenotypes = Vector::Zero( 1 ) ;
	
		snptest2::NullModelLogLikelihood ll( phenotypes ) ;
		Vector parameters = Vector::Zero( 1 ) ;
		ll.evaluate_at( parameters ) ;
		TEST_ASSERT( ll.get_value_of_function() == std::log( 0.5 ) ) ;
		TEST_ASSERT( ll.get_value_of_first_derivative() == Vector::Constant( 1, -0.5 ) ) ;
		TEST_ASSERT( ll.get_value_of_second_derivative() == Matrix::Constant( 1, 1, -0.25 )) ;

		parameters( 0 ) = 1.0 ;
		double p0 = 1.0 / ( 1.0 + std::exp( 1.0 )) ;
		double p1 = std::exp(1.0) / ( 1.0 + std::exp( 1.0 )) ;
		ll.evaluate_at( parameters ) ;
		TEST_ASSERT( ll.get_value_of_function() == std::log( p0 )) ;
		TEST_ASSERT( ll.get_value_of_first_derivative() == Vector::Constant( 1, -p1 )) ;
		TEST_ASSERT( ll.get_value_of_second_derivative() == Matrix::Constant( 1, 1, -p1 * p0 )) ;
	}
	
	{
		Vector phenotypes( 2 ) ;
		phenotypes << 0.0, 1.0 ;
		snptest2::NullModelLogLikelihood ll( phenotypes ) ;
		Vector parameters = Vector::Zero( 1 ) ;
		ll.evaluate_at( parameters ) ;
		TEST_ASSERT( ll.get_value_of_function() == std::log( 0.5 ) + std::log( 0.5 ) ) ;
		TEST_ASSERT( ll.get_value_of_first_derivative() == Vector::Constant( 1, 0.0 ) ) ;
		TEST_ASSERT( ll.get_value_of_second_derivative() == Matrix::Constant( 1, 1, -0.5 )) ;
		
		parameters( 0 ) = 1.0 ;
		double p0 = 1.0 / ( 1.0 + std::exp( 1.0 )) ;
		double p1 = std::exp(1.0) / ( 1.0 + std::exp( 1.0 )) ;
		ll.evaluate_at( parameters ) ;
		TEST_ASSERT( ll.get_value_of_function() == std::log( p0 ) + std::log( p1 )) ;
		TEST_ASSERT( ll.get_value_of_first_derivative() == Vector::Constant( 1, 1.0 - 2.0 * p1 )) ;
		TEST_ASSERT( ll.get_value_of_second_derivative() == Matrix::Constant( 1, 1, -2.0 * p1 * ( 1.0 - p1 ) )) ;
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_MAIN
{
	test_null_model_small_datasets() ;
}
