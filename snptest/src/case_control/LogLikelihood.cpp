#include <iostream>
#include <vector>
#include <utility>
#include <boost/iterator/counting_iterator.hpp>
#include "snptest/case_control/LogLikelihood.hpp"

namespace snptest {
	namespace case_control {
		LogLikelihood::LogLikelihood(
			Vector const& phenotypes,
			FinitelySupportedFunctionSet const& genotypes
		):
			m_phenotypes( phenotypes ),
			m_genotypes( genotypes ),
			m_design_matrices( calculate_design_matrices( genotypes.get_support() ) ),
			m_included_samples(
				std::vector< std::size_t >(
					boost::counting_iterator< std::size_t >( 0 ),
					boost::counting_iterator< std::size_t >( m_phenotypes.size() )
				)
			)
		{
			assert( m_phenotypes.size() == m_genotypes.get_number_of_functions() ) ;
			assert( m_included_samples.size() <= std::size_t( m_genotypes.get_number_of_functions() ) ) ;
			assert( m_design_matrices.size() == 3 ) ;
		}

		LogLikelihood::LogLikelihood(
			Vector const& phenotypes,
			FinitelySupportedFunctionSet const& genotypes,
			std::vector< std::size_t > const& included_samples
		):
			m_phenotypes( phenotypes ),
			m_genotypes( genotypes ),
			m_design_matrices( calculate_design_matrices( genotypes.get_support() ) ),
			m_included_samples( included_samples )
		{
			assert( m_phenotypes.size() == m_genotypes.get_number_of_functions() ) ;
			assert( m_included_samples.size() <= std::size_t( m_genotypes.get_number_of_functions() ) ) ;
			assert( m_design_matrices.size() == 3 ) ;
		}

		std::vector< LogLikelihood::Vector > LogLikelihood::calculate_design_matrices(
			Vector const& genotype_levels
		) const {
			assert( genotype_levels.size() == 3 ) ;
			std::vector< Vector > result( 3, Vector::Zero( 2 ) ) ;
			for( std::size_t g = 0; g < 3; ++g ) {
				result[g]( 0 ) = 1.0 ;
				result[g]( 1 ) = genotype_levels(g) ;
			}
			return result ;
		}

		// Calculate  p_{g,\theta}(\psi_i) for each individual i and each genotype, given the parameters.
		// (This is the probability of the given phenotype from logistic regression).
		LogLikelihood::Matrix LogLikelihood::calculate_p_g_thetas(
			Vector const& parameters,
			std::vector< Vector > const& design_matrices
		) const {
			assert( design_matrices.size() == 3 ) ;
			Matrix result( 2, 3 ) ;
			for( std::size_t g = 0; g < 3; ++g ) {
				result( 0, g ) = calculate_p_g_theta( parameters, design_matrices[g]( 1 ), 0.0 ) ;
				result( 1, g ) = calculate_p_g_theta( parameters, design_matrices[g]( 1 ), 1.0 ) ;
				assert( result( 0, g ) == result( 0, g ) ) ;
				assert( result( 1, g ) == result( 1, g ) ) ;
			}
			return result ;
		}

		void LogLikelihood::evaluate_at( Point const& parameters ) {
			assert( m_design_matrices.size() == 3 ) ;
			m_p_g_thetas = calculate_p_g_thetas( parameters, m_design_matrices ) ;
			m_coefficients = calculate_coefficients( m_genotypes.get_values(), m_phenotypes, m_p_g_thetas ) ;
		}
	
		// Calculate the z coefficients, i.e.
		// z_ig = p( g_i = g | \psi_i, o_i, \theta )
		// By Bayes theorem this is given as
		//              p_{g,\theta}( \psi_i ) c_{ig}
		//   z_ig =  ---------------------------------
		//            \sum_g p_g,\theta( \psi_i ) c_ig
		//
		// where c_ig is the call probability.
		LogLikelihood::Matrix LogLikelihood::calculate_coefficients(
			Matrix const& genotypes,
			Vector const& phenotypes,
			Matrix const& p_g_thetas
		) const {
			Matrix result = genotypes ;

			for( int i = 0; i < result.rows(); ++i ) {
				result.row( i ) = result.row( i ).cwiseProduct( p_g_thetas.row( std::size_t( phenotypes( i )) )) ;
				double sum = result.row( i ).sum() ;
				if( sum > 0.0 ) {
					result.row( i ) = result.row( i ) / result.row( i ).sum() ;
				}
				else {
					result.row( i ) = Vector::Constant( 3, std::numeric_limits< double >::quiet_NaN() ) ;
				}
			}

			for( int i = 0; i < phenotypes.size(); ++i ) {
				for( std::size_t g = 0; g < 3; ++g ) {
					result( i, g ) *= ( phenotypes( i ) - p_g_thetas( 1, g ) ) ;
				}
			}

			return result ;
		}

		double LogLikelihood::get_value_of_function() const {
			return calculate_function(
				m_phenotypes,
				m_genotypes.get_values(),
				m_p_g_thetas
			) ;
		}

		LogLikelihood::Vector LogLikelihood::get_value_of_first_derivative() const {
			return calculate_first_derivative(
				m_phenotypes,
				m_p_g_thetas,
				m_coefficients,
				m_design_matrices
			) ;
		}

		LogLikelihood::Matrix LogLikelihood::get_value_of_second_derivative() const {
			return calculate_second_derivative(
				m_phenotypes,
				m_p_g_thetas,
				m_coefficients,
				m_design_matrices
			) ;
		}

		double LogLikelihood::calculate_function(
			Vector const& phenotypes,
			Matrix const& genotypes,
			Matrix const& p_g_thetas
		) const {
			double result = 0.0 ;
			for( std::size_t i = 0; i < m_included_samples.size(); ++i ) {
				int const sample_i = m_included_samples[i] ;
				double x = 0.0 ;
				for( std::size_t g = 0; g < 3; ++g ) {
					x += genotypes( sample_i, g ) * p_g_thetas( int( phenotypes( sample_i )), g ) ;
				}
				result += std::log( x ) ;
			}
			return result ;
		}

		LogLikelihood::Vector LogLikelihood::calculate_first_derivative(
			Vector const& phenotypes,
			Matrix const& p_g_thetas,
			Matrix const& coefficients,
			std::vector< Vector > const& design_matrices
		) const {
			Vector result = Vector::Zero( 2 ) ;
			for( std::size_t i = 0; i < m_included_samples.size(); ++i ) {
				int const sample_i = m_included_samples[i] ;
				for( std::size_t g = 0; g < 3; ++g ) {
					result += coefficients( sample_i, g ) * design_matrices[g] ;
				}
			}
			return result ;
		}

		LogLikelihood::Matrix LogLikelihood::calculate_second_derivative(
			Vector const& phenotypes,
			Matrix const& p_g_thetas,
			Matrix const& coefficients,
			std::vector< Vector > const& design_matrices
		) const {
			Matrix result = Matrix::Zero( 2, 2 ) ;
			for( std::size_t i = 0; i < m_included_samples.size(); ++i ) {
				int const sample_i = m_included_samples[i] ;
				Vector second_term = Vector::Zero( 2 ) ;
				for( std::size_t g = 0; g < 3; ++g ) {
					result += coefficients( sample_i, g ) * ( 1.0 - 2.0 * p_g_thetas( 1, g )) * ( design_matrices[g] * design_matrices[g].transpose() ) ;
					second_term += coefficients( sample_i, g ) * design_matrices[g] ;
				}
				result -= second_term * second_term.transpose() ;
			}
			return result ;
		}
		
	}
}
