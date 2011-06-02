#include <iostream>
#include <iomanip>
#include "test_case.hpp"
#include "snptest/SNPTEST2AlternativeModel.hpp"
#include "Eigen/Eigen"

void test_alternative_model_certain_genotypes_one_individual( std::size_t g ) {
	std::cerr << "Testing SNPTEST2 alternative model (one individual with certain genotype " << g << ")..." ;
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::MatrixXd Matrix ;

	std::vector< double > genotype_levels( 3 ) ;
	genotype_levels[0] = 0.0 ;
	genotype_levels[1] = 1.0 ;
	genotype_levels[2] = 2.0 ;

	Vector const phenotypes = Vector::Zero( 1 ) ;
	Matrix genotypes = Matrix::Zero( 1, 3 )  ;
	genotypes( 0, g ) = 1.0 ;

	Vector design_matrix( 2 ) ;
	design_matrix << 1.0, g ;

	Matrix design_matrix_squared( 2, 2 ) ;
	design_matrix_squared <<	1.0, g,
								g, g*g ;

	double p0 ;
	double p1 ;

	snptest2::AlternativeModelLogLikelihood ll( phenotypes, genotypes, genotype_levels ) ;
	Vector parameters( 2 ) ;
	parameters << 0.0, 0.0 ;

	p1 = std::exp( parameters(0) + genotype_levels[g] * parameters(1) ) ;
	p0 = 1.0 / ( 1.0 + p1 ) ;
	p1 = p1 / ( 1.0 + p1 ) ;

	ll.evaluate_at( parameters ) ;
	TEST_ASSERT( ll.get_value_of_function() == std::log( p0 ) ) ;
	TEST_ASSERT( ll.get_value_of_first_derivative() == -p1 * design_matrix ) ;
	TEST_ASSERT( ll.get_value_of_second_derivative() == (( -p1 * ( 1.0 - 2.0 * p1 )) - ( p1 * p1 )) * design_matrix_squared ) ;

	parameters << 0.0, 1.0 ;
	p1 = std::exp( parameters(0) + genotype_levels[g] * parameters(1) ) ;
	p0 = 1.0 / ( 1.0 + p1 ) ;
	p1 = p1 / ( 1.0 + p1 ) ;
	
	ll.evaluate_at( parameters ) ;
	TEST_ASSERT( ll.get_value_of_function() == std::log( p0 ) ) ;
	TEST_ASSERT( ll.get_value_of_first_derivative() == -p1 * design_matrix ) ;
	TEST_ASSERT( ll.get_value_of_second_derivative() == (( -p1 * ( 1.0 - 2.0 * p1 )) - ( p1 * p1 )) * design_matrix_squared ) ;

	// Change baseline parameter
	parameters << 1.0, 0.0 ;
	p1 = std::exp( parameters(0) + genotype_levels[g] * parameters(1) ) ;
	p0 = 1.0 / ( 1.0 + p1 ) ;
	p1 = p1 / ( 1.0 + p1 ) ;
	ll.evaluate_at( parameters ) ;
	TEST_ASSERT( ll.get_value_of_function() == std::log( p0 ) ) ;
	TEST_ASSERT( ll.get_value_of_first_derivative() == -p1 * design_matrix ) ;
	TEST_ASSERT( ll.get_value_of_second_derivative() == (( -p1 * ( 1.0 - 2.0 * p1 )) - ( p1 * p1 )) * design_matrix_squared ) ;

	parameters << 1.0, 1.0 ;
	p1 = std::exp( parameters(0) + genotype_levels[g] * parameters(1) ) ;
	p0 = 1.0 / ( 1.0 + p1 ) ;
	p1 = p1 / ( 1.0 + p1 ) ;
	
	ll.evaluate_at( parameters ) ;
	TEST_ASSERT( ll.get_value_of_function() == std::log( p0 ) ) ;
	TEST_ASSERT( ll.get_value_of_first_derivative() == -p1 * design_matrix ) ;
	TEST_ASSERT( ll.get_value_of_second_derivative() == (( -p1 * ( 1.0 - 2.0 * p1 )) - ( p1 * p1 )) * design_matrix_squared ) ;

	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_alternative_model_one_individual ) {
	test_alternative_model_certain_genotypes_one_individual( 0 ) ;
	test_alternative_model_certain_genotypes_one_individual( 1 ) ;
	test_alternative_model_certain_genotypes_one_individual( 2 ) ;
}

