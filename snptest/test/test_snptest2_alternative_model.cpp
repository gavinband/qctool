
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include "test_case.hpp"
#include "snptest/case_control/AlternativeModelLogLikelihood.hpp"
#include "Eigen/Eigen"

typedef Eigen::VectorXd Vector ;
typedef Eigen::MatrixXd Matrix ;

namespace {
	template< typename M1, typename M2 >
	void check_equal( M1 const& left, M2 const& right, double tolerance = 0.0001 ) {
		BOOST_CHECK_EQUAL( left.rows(), right.rows() ) ;
		BOOST_CHECK_EQUAL( left.cols(), right.cols() ) ;
		for( int i = 0; i < left.rows(); ++i ) {
			for( int j = 0; j < left.cols(); ++j ) {
				BOOST_CHECK_CLOSE( left(i,j), right(i,j), tolerance ) ;
			}
		}
	}
}

void test_alternative_model_certain_genotypes_one_individual( std::size_t g ) {
	std::cerr << "Testing SNPTEST2 alternative model (one individual with certain genotype " << g << ")..." ;
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::MatrixXd Matrix ;

	Vector const phenotypes = Vector::Zero( 1 ) ;
	Matrix genotypes = Matrix::Zero( 1, 3 )  ;
	Vector genotype_levels( 3 ) ;
	genotype_levels << 0, 1, 2 ;
	genotypes( 0, g ) = 1.0 ;

	Vector design_matrix( 2 ) ;
	design_matrix << 1.0, g ;

	Matrix design_matrix_squared( 2, 2 ) ;
	design_matrix_squared <<	1.0, g,
								g, g*g ;

	double p0 ;
	double p1 ;

	snptest::case_control::AlternativeModelLogLikelihood ll ;
	ll.set_phenotypes( phenotypes ) ;
	ll.set_predictor_probs( genotypes, genotype_levels ) ;

	Vector parameters( 2 ) ;
	// Start at 0,0
	parameters << 0.0, 0.0 ;

	p1 = std::exp( parameters(0) + genotype_levels(g) * parameters(1) ) ;
	p0 = 1.0 / ( 1.0 + p1 ) ;
	p1 = p1 / ( 1.0 + p1 ) ;

	ll.evaluate_at( parameters ) ;
	BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( p0 ) ) ;
	BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), -p1 * design_matrix ) ;
	BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), (( -p1 * ( 1.0 - 2.0 * p1 )) - ( p1 * p1 )) * design_matrix_squared ) ;

	// Nonzero genetic effect parameter
	parameters << 0.0, 1.0 ;
	p1 = std::exp( parameters(0) + genotype_levels(g) * parameters(1) ) ;
	p0 = 1.0 / ( 1.0 + p1 ) ;
	p1 = p1 / ( 1.0 + p1 ) ;
	
	ll.evaluate_at( parameters ) ;
	BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( p0 ) ) ;
	check_equal( ll.get_value_of_first_derivative(), -p1 * design_matrix, 0.00000000001 ) ;
	BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), (( -p1 * ( 1.0 - 2.0 * p1 )) - ( p1 * p1 )) * design_matrix_squared ) ;

	// Nonzero baseline parameter instead
	parameters << 1.0, 0.0 ;
	p1 = std::exp( parameters(0) + genotype_levels(g) * parameters(1) ) ;
	p0 = 1.0 / ( 1.0 + p1 ) ;
	p1 = p1 / ( 1.0 + p1 ) ;
	ll.evaluate_at( parameters ) ;
	BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( p0 ) ) ;
	BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), -p1 * design_matrix ) ;
	BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), (( -p1 * ( 1.0 - 2.0 * p1 )) - ( p1 * p1 )) * design_matrix_squared ) ;

	// Both parameters nonzero
	parameters << 1.0, 1.0 ;
	p1 = std::exp( parameters(0) + genotype_levels(g) * parameters(1) ) ;
	p0 = 1.0 / ( 1.0 + p1 ) ;
	p1 = p1 / ( 1.0 + p1 ) ;
	
	ll.evaluate_at( parameters ) ;
	BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( p0 ) ) ;
	check_equal( ll.get_value_of_first_derivative(), -p1 * design_matrix, 0.000000000001 ) ;
	Matrix expected_2nd_derivative = (( -p1 * ( 1.0 - 2.0 * p1 )) - ( p1 * p1 )) * design_matrix_squared ;
	BOOST_CHECK_SMALL(( ll.get_value_of_second_derivative() - expected_2nd_derivative ).maxCoeff(), 0.0000000001 ) ;

	std::cerr << "ok.\n" ;
}

namespace {
	double mean_fn( Vector const parameters, double phenotype, double predictor ) {
		double result = std::exp( parameters(0) + predictor * parameters( 1 )) ;
		if( phenotype == 1 ) {
			return result / ( 1.0 + result ) ;
		}
		else {
			return 1.0 / ( 1.0 + result ) ;
		}
	}

	void tamugoi(
		Vector phenotypes,
		Matrix const& genotypes,
		Vector const& genotype_levels,
		Vector const& parameters,
		std::size_t g1,
		std::size_t g2
	) {
		snptest::case_control::AlternativeModelLogLikelihood ll ;
		ll.set_phenotypes( phenotypes ) ;
		ll.set_predictor_probs( genotypes, genotype_levels ) ;

		Matrix design_matrix1( 1, 2 ) ;
		design_matrix1.col(0).setOnes() ;
		design_matrix1( 0, 1 ) = genotype_levels( g1 ) ;
		Matrix design_matrix2( 1, 2 ) ;
		design_matrix2.col(0).setOnes() ;
		design_matrix2( 0, 1 ) = genotype_levels( g2 ) ;
	
		double const f1 = mean_fn( parameters, phenotypes(0), genotype_levels( g1 ) ) ;
		double const f2 = mean_fn( parameters, phenotypes(0), genotype_levels( g2 ) ) ;

		ll.evaluate_at( parameters ) ;
		double h1
			= genotypes( 0, g1 ) * f1
			+ genotypes( 0, g2 ) * f2 ;
	
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), std::log( h1 ) ) ;

		Vector Dh1
			= ( genotypes( 0, g1 ) * f1 * ( 1 - f1 ) * design_matrix1.row(0).transpose() )
			+ ( genotypes( 0, g2 ) * f2 * ( 1 - f2 ) * design_matrix2.row(0).transpose() ) ;
	
		if( phenotypes(0) == 0.0 ) {
			Dh1 = -Dh1 ;
		}
	
		BOOST_CHECK_SMALL( ( ll.get_value_of_first_derivative() - ( Dh1/h1 ) ).array().abs().maxCoeff(), 0.0000000001 ) ;

		Matrix D2h1
			= ( genotypes( 0, g1 ) * f1 * ( 1 - f1 ) * ( 1 - 2*f1 ) * design_matrix1.row(0).transpose() * design_matrix1.row(0) )
			+ ( genotypes( 0, g2 ) * f2 * ( 1 - f2 ) * ( 1 - 2*f2 ) * design_matrix2.row(0).transpose() * design_matrix2.row(0) ) ;
		
		Matrix adjustment = (Dh1/h1) * (Dh1/h1).transpose() ;
		BOOST_CHECK_SMALL( ( ll.get_value_of_second_derivative() - ((D2h1/h1) - adjustment) ).array().abs().maxCoeff(), 0.0000000001 ) ;
	}
}

void test_alternative_model_uncertain_genotypes_one_individual( std::size_t g1, std::size_t g2, double v ) {
	std::cerr << "test_alternative_model_uncertain_genotypes_one_individual with genotypes " << g1 << "(" << v << "), " << g2 << "(" << (1.0-v)  << ")..." ;

	Vector const phenotypes = Vector::Zero( 1 ) ;
	Matrix genotypes = Matrix::Zero( 1, 3 )  ;
	Vector genotype_levels( 3 ) ;
	genotype_levels << 0, 1, 2 ;
	genotypes( 0, g1 ) = v ;
	genotypes( 0, g2 ) += ( 1.0 - v ) ;

	Vector parameters( 2 ) ;
	// Start at 0,0
	parameters << 0.0, 0.0 ;
	tamugoi( phenotypes, genotypes, genotype_levels, parameters, g1, g2 ) ;
	
	parameters << 0.0, 1.0 ;
	tamugoi( phenotypes, genotypes, genotype_levels, parameters, g1, g2 ) ;

	parameters << 1.0, 0.0 ;
	tamugoi( phenotypes, genotypes, genotype_levels, parameters, g1, g2 ) ;

	parameters << 1.0, 1.0 ;
	tamugoi( phenotypes, genotypes, genotype_levels, parameters, g1, g2 ) ;

//	parameters << 50.0, 24.0 ;
//	tamugoi( phenotypes, genotypes, genotype_levels, parameters, g1, g2 ) ;
	
	double mean_genotype = g1 * genotypes( 0, g1 ) + g2 * genotypes( 0, g2 ) ;
	genotype_levels << -mean_genotype, 1.0 - mean_genotype, 2.0 - mean_genotype ;

	parameters << 0.0, 0.0 ;
	tamugoi( phenotypes, genotypes, genotype_levels, parameters, g1, g2 ) ;
	
	parameters << 0.0, 1.0 ;
	tamugoi( phenotypes, genotypes, genotype_levels, parameters, g1, g2 ) ;

	parameters << 1.0, 0.0 ;
	tamugoi( phenotypes, genotypes, genotype_levels, parameters, g1, g2 ) ;

	parameters << 1.0, 1.0 ;
	tamugoi( phenotypes, genotypes, genotype_levels, parameters, g1, g2 ) ;

//	parameters << 50.0, 24.0 ;
//	tamugoi( phenotypes, genotypes, genotype_levels, parameters, g1, g2 ) ;
}


AUTO_TEST_CASE( test_alternative_model_one_individual ) {
	test_alternative_model_certain_genotypes_one_individual( 0 ) ;
	test_alternative_model_certain_genotypes_one_individual( 1 ) ;
	test_alternative_model_certain_genotypes_one_individual( 2 ) ;
	
	test_alternative_model_uncertain_genotypes_one_individual( 0, 1, 0.5 ) ;
	test_alternative_model_uncertain_genotypes_one_individual( 0, 2, 0.5 ) ;
	test_alternative_model_uncertain_genotypes_one_individual( 1, 2, 0.5 ) ;
	test_alternative_model_uncertain_genotypes_one_individual( 0, 1, 0.2 ) ;
	test_alternative_model_uncertain_genotypes_one_individual( 0, 2, 0.2 ) ;
	test_alternative_model_uncertain_genotypes_one_individual( 1, 2, 0.2 ) ;
}

