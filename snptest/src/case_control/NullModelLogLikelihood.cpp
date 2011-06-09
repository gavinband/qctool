#include <iostream>
#include <vector>
#include <boost/iterator/counting_iterator.hpp>
#include "snptest/case_control/NullModelLogLikelihood.hpp"

namespace snptest {
	namespace case_control {
		NullModelLogLikelihood::NullModelLogLikelihood(
			Vector const& phenotypes,
			FinitelySupportedFunctionSet const& genotypes,
			bool weight_by_genotypes
		):
			m_phenotypes( phenotypes ),
			m_genotypes( genotypes ),
			m_weight_by_genotypes( weight_by_genotypes ),
			m_included_samples( std::vector< std::size_t >(
					boost::counting_iterator< std::size_t >( 0 ),
					boost::counting_iterator< std::size_t >( phenotypes.size() )
				)
			)
		{}

		NullModelLogLikelihood::NullModelLogLikelihood(
			Vector const& phenotypes,
			FinitelySupportedFunctionSet const& genotypes,
			bool weight_by_genotypes,
			std::vector< std::size_t > const& included_samples
		):
			m_phenotypes( phenotypes ),
			m_genotypes( genotypes ),
			m_weight_by_genotypes( weight_by_genotypes ),
			m_included_samples( included_samples )
		{
			assert( m_included_samples.size() <= std::size_t( m_phenotypes.size() ) ) ;
		}

		void NullModelLogLikelihood::evaluate_at( Vector const& parameters ) {
			m_p_thetas = calculate_p_thetas( parameters ) ;
		}

		// Calculate  p_{\theta}(\psi_i) for each individual i, given the parameters.
		// (This is the probability of the given phenotype from logistic regression).
		NullModelLogLikelihood::Matrix NullModelLogLikelihood::calculate_p_thetas(
			Vector const& parameters
		) const {
			assert( parameters.size() == 1 ) ;
			Matrix result( 2, 1 ) ;
			result( 0, 0 ) = calculate_p_theta( parameters, 0.0 ) ;
			result( 1, 0 ) = calculate_p_theta( parameters, 1.0 ) ;
			return result ;
		}

		double NullModelLogLikelihood::calculate_p_theta(
			Vector const& parameters,
			double phenotype
		) const {
			double x = std::exp( parameters(0) ) ;
			return (( phenotype == 1.0 ) ? x : 1.0 ) / ( 1 + x ) ;
		}
	
		double NullModelLogLikelihood::get_value_of_function() const {
			return calculate_function(
				m_phenotypes,
				m_p_thetas
			) ;
		}

		double NullModelLogLikelihood::calculate_function(
			Vector const& phenotypes,
			Matrix const& p_thetas
		) const {
			double result = 0.0 ;
			for( std::size_t i = 0; i < m_included_samples.size(); ++i ) {
				int const sample_i = m_included_samples[i] ;
				result += std::log(
					( m_weight_by_genotypes ? m_genotypes.get_values( sample_i ).sum() : 1.0 ) * p_thetas( int( phenotypes( sample_i ) ) ) ) ;
			}
			return result ;
		}

		NullModelLogLikelihood::Vector NullModelLogLikelihood::get_value_of_first_derivative() const {
			return calculate_first_derivative(
				m_phenotypes,
				m_p_thetas
			) ;
		}
	
		NullModelLogLikelihood::Matrix NullModelLogLikelihood::get_value_of_second_derivative() const {
			return calculate_second_derivative(
				m_phenotypes,
				m_p_thetas
			) ;
		}
	
		NullModelLogLikelihood::Vector NullModelLogLikelihood::calculate_first_derivative(
			Vector const& phenotypes,
			Matrix const& p_thetas
		) const {
			Vector result = Vector::Zero( 1 ) ;
			for( std::size_t i = 0; i < m_included_samples.size(); ++i ) {
				int const sample_i = m_included_samples[i] ;
				result( 0 ) += phenotypes( sample_i ) - p_thetas( 1 ) ;
				//result( 0 ) = result(0) + phenotypes( i ) - p_thetas( 1 ) ;
			}
			return result ;
		}
	
		NullModelLogLikelihood::Matrix NullModelLogLikelihood::calculate_second_derivative(
			Vector const& phenotypes,
			Matrix const& p_thetas
		) const {
			return Matrix::Constant(
				1, 1,
				-double( phenotypes.size() ) * p_thetas( 1 ) * ( 1.0 - p_thetas( 1 ))
			) ;
		}
	}
}
