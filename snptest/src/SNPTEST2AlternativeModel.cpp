#include <iostream>
#include <vector>
#include <utility>
#include "snptest/SNPTEST2AlternativeModel.hpp"

namespace snptest2 {
	AlternativeModelLogLikelihood::AlternativeModelLogLikelihood(
		AlternativeModelLogLikelihood const& other
	):
		m_phenotypes( other.m_phenotypes ),
		m_genotypes( other.m_genotypes ),
		m_design_matrices( other.m_design_matrices ),
		m_p_g_thetas( other.m_p_g_thetas ),
		m_coefficients( other.m_coefficients )
	{}

	AlternativeModelLogLikelihood::AlternativeModelLogLikelihood(
		Vector const& phenotypes,
		Matrix const& genotypes,
		std::vector< double > const& genotype_levels
	):
		m_phenotypes( phenotypes ),
		m_genotypes( genotypes ),
		m_design_matrices( calculate_design_matrices( genotype_levels ) )
	{
		assert( m_phenotypes.size() == m_genotypes.rows() ) ;
		assert( m_design_matrices.size() == 3 ) ;
	}

	std::vector< AlternativeModelLogLikelihood::Vector > AlternativeModelLogLikelihood::calculate_design_matrices(
		std::vector< double > const& genotype_levels
	) const {
		assert( genotype_levels.size() == 3 ) ;
		std::vector< Vector > result( 3, Vector::Zero( 2 ) ) ;
		for( std::size_t g = 0; g < 3; ++g ) {
			result[g]( 0 ) = 1.0 ;
			result[g]( 1 ) = genotype_levels[g] ;
		}
		return result ;
	}

	// Calculate  p_{g,\theta}(\psi_i) for each individual i and each genotype, given the parameters.
	// (This is the probability of the given phenotype from logistic regression).
	AlternativeModelLogLikelihood::Matrix AlternativeModelLogLikelihood::calculate_p_g_thetas(
		Vector const& parameters,
		std::vector< Vector > const& design_matrices
	) const {
		assert( parameters.size() == 2 ) ;
		assert( design_matrices.size() == 3 ) ;
		Matrix result( 2, 3 ) ;
		for( std::size_t g = 0; g < 3; ++g ) {
			result( 0, g ) = calculate_p_g_theta( parameters, design_matrices[g]( 1 ), 0.0 ) ;
			result( 1, g ) = calculate_p_g_theta( parameters, design_matrices[g]( 1 ), 1.0 ) ;
		}
		return result ;
	}

	double AlternativeModelLogLikelihood::calculate_p_g_theta(
		Vector const& parameters,
		double genotype,
		double phenotype
	) const {
		double x = std::exp( parameters(0) + parameters(1) * genotype ) ;
		return (( phenotype == 1.0 ) ? x : 1.0 ) / ( 1 + x ) ;
	}

	void AlternativeModelLogLikelihood::evaluate_at( Point const& parameters ) {
		assert( parameters.size() == 2 ) ;
		assert( m_design_matrices.size() == 3 ) ;
		m_p_g_thetas = calculate_p_g_thetas( parameters, m_design_matrices ) ;
		m_coefficients = calculate_coefficients( m_genotypes, m_phenotypes, m_p_g_thetas ) ;
	}
	
	double AlternativeModelLogLikelihood::get_value_of_function() const {
		return calculate_function(
			m_phenotypes,
			m_genotypes,
			m_p_g_thetas
		) ;
	}

	AlternativeModelLogLikelihood::Vector AlternativeModelLogLikelihood::get_value_of_first_derivative() const {
		return calculate_first_derivative(
			m_phenotypes,
			m_p_g_thetas,
			m_coefficients,
			m_design_matrices
		) ;
	}

	AlternativeModelLogLikelihood::Matrix AlternativeModelLogLikelihood::get_value_of_second_derivative() const {
		return calculate_second_derivative(
			m_phenotypes,
			m_p_g_thetas,
			m_coefficients,
			m_design_matrices
		) ;
	}

	// Calculate the z coefficients, i.e.
	// z_ig = p( g_i = g | \psi_i, o_i, \theta )
	// By Bayes theorem this is given as
	//              p_{g,\theta}( \psi_i ) c_{ig}
	//   z_ig =  ---------------------------------
	//            \sum_g p_g,\theta( \psi_i ) c_ig
	//
	// where c_ig is the call probability.
	AlternativeModelLogLikelihood::Matrix AlternativeModelLogLikelihood::calculate_coefficients(
		Matrix const& genotypes,
		Vector const& phenotypes,
		Matrix const& p_g_thetas
	) const {
		Matrix result = genotypes ;

		for( int i = 0; i < result.rows(); ++i ) {
			result.row( i ) = result.row( i ).cwiseProduct( p_g_thetas.row( std::size_t( phenotypes( i )) )) ;
			result.row( i ) = result.row( i ) / result.row( i ).sum() ;
		}

		for( int i = 0; i < phenotypes.size(); ++i ) {
			for( std::size_t g = 0; g < 3; ++g ) {
				result( i, g ) *= ( phenotypes( i ) - p_g_thetas( 1, g ) ) ;
			}
		}
		return result ;
	}

	double AlternativeModelLogLikelihood::calculate_function(
		Vector const& phenotypes,
		Matrix const& m_genotypes,
		Matrix const& p_g_thetas
	) const {
		double result = 0.0 ;
		for( int i = 0; i < phenotypes.size(); ++i ) {
			double x = 0.0 ;
			for( std::size_t g = 0; g < 3; ++g ) {
				x += m_genotypes( i, g ) * p_g_thetas( int( phenotypes( i )), g ) ;
			}
			result += std::log( x ) ;
		}
		return result ;
	}

	AlternativeModelLogLikelihood::Vector AlternativeModelLogLikelihood::calculate_first_derivative(
		Vector const& phenotypes,
		Matrix const& p_g_thetas,
		Matrix const& coefficients,
		std::vector< Vector > const& design_matrices
	) const {
		Vector result = Vector::Zero( 2 ) ;
		for( int i = 0; i < phenotypes.size(); ++i ) {
			for( std::size_t g = 0; g < 3; ++g ) {
				result += coefficients( i, g ) * design_matrices[g] ;
			}
		}
		return result ;
	}

	AlternativeModelLogLikelihood::Matrix AlternativeModelLogLikelihood::calculate_second_derivative(
		Vector const& phenotypes,
		Matrix const& p_g_thetas,
		Matrix const& coefficients,
		std::vector< Vector > const& design_matrices
	) const {
		Matrix result = Matrix::Zero( 2, 2 ) ;
		for( int i = 0; i < phenotypes.size(); ++i ) {
			Vector second_term = Vector::Zero( 2 ) ;
			for( std::size_t g = 0; g < 3; ++g ) {
				result += coefficients(i, g) * ( 1.0 - 2.0 * p_g_thetas( 1, g )) * ( design_matrices[g] * design_matrices[g].transpose() ) ;
				second_term += coefficients(i, g) * design_matrices[g] ;
			}
			result -= second_term * second_term.transpose() ;
		}
		return result ;
	}
	
}