#include <iostream>
#include <vector>
#include "snptest/SNPTEST2NullModel.hpp"

namespace snptest2 {
	NullModelLogLikelihood::NullModelLogLikelihood(
		NullModelLogLikelihood const& other
	):
		m_phenotypes( other.m_phenotypes )
	{
		assert(0) ;
	}

	NullModelLogLikelihood::NullModelLogLikelihood(
		Vector const& phenotypes
	):
		m_phenotypes( phenotypes )
	{}

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
		for( int i = 0; i < phenotypes.rows(); ++i ) {
			result += std::log( p_thetas( int( phenotypes(i) ) ) ) ;
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
		for( int i = 0; i < phenotypes.size(); ++i ) {
			result( 0 ) += phenotypes( i ) - p_thetas( 1 ) ;
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
