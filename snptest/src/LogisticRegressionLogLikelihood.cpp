
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <vector>
#include <utility>
#include "snptest/LogisticRegressionLogLikelihood.hpp"
#include "snptest/FinitelySupportedFunctionSet.hpp"

namespace snptest2 {
	LogisticRegressionLogLikelihood::LogisticRegressionLogLikelihood(
		LogisticRegressionLogLikelihood const& other
	):
		m_outcomes( other.m_outcomes ),
		m_design_matrices( other.m_design_matrices ),
		m_design_probabilities( other.m_design_probabilities ),
		m_p_outcome_given_predictors_and_parameters( other.m_p_outcome_given_predictors_and_parameters ),
		m_coefficients( other.m_coefficients )
	{}

	LogisticRegressionLogLikelihood::LogisticRegressionLogLikelihood(
		Vector const& outcomes,
		Matrix const& design_matrices,
		Matrix const& design_probabilities
	):
		m_outcomes( outcomes ),
		m_design_matrices( design_matrices ),
		m_design_probabilities( design_probabilities )
	{
		assert( m_design_probabilities.rows() == outcomes.size() ) ;
		assert( m_design_matrices.cols() == m_design_probabilities.cols() ) ;
	}

	#if 0
	void LogisticRegressionLogLikelihood::calculate_design(
		std::vector< FinitelySupportedFunctionSet > const& predictors,
		Matrix* design_matrices,
		Matrix* design_probabilities
	) const {
		std::size_t const N = predictors.size() ;

		// Get number of predictor levels.
		std::size_t ncols = predictors.front().get_size_of_support() ;
		for( std::size_t i = 1; i < predictors.size(); ++i ) {
			ncols *= predictors[i].get_size_of_support() ;
		}
		
		design_matrices->resize( predictors.size() + 1, ncols ) ;
		design_matrices->row(0) = Vector::Ones( ncols ) ;
		design_probabilities->resize( predictors.size() + 1, ncols ) ;

		std::vector< std::size_t > choice( N, 0ul ) ;
		bool finished = ( ncols == 0 ) ;
		std::size_t count = 0 ;
		for( ; !finished; ++count ) {
			// Store the design matrix and predictor probabilities for this choice of predictors.
			design_matrices->operator()( 1, count ) = predictors[0].get_support()( choice[0] ) ;
			design_probabilities->col( count ) = predictors[0].get_values().col( choice[ 0 ] ) ;
			for( std::size_t i = 1; i < predictors.size(); ++i ) {
				design_matrices->operator()( i+1, count ) = predictors[i].get_support()( choice[i] ) ;
				design_probabilities->col( count ).array() *= predictors[i].get_values().col( choice[ i ] ).array() ;
			}

			// look at next possible combination.
			std::size_t current_elt = N - 1 ;
			for( ; ++choice[ current_elt ] == predictors[ current_elt ].get_size_of_support() ; --current_elt ) {
				if( current_elt > 0 ) {
					choice[ current_elt ] = 0 ;
				}
				else {
					finished = true ;
				}
			}
		}
		assert( count == ncols ) ;
	}
	#endif
	
	void LogisticRegressionLogLikelihood::evaluate_at( Point const& parameters ) {
		assert( parameters.size() == m_design_matrices.rows() ) ;
		calculate_p_outcome_given_predictors_and_parameters( parameters, m_design_matrices, &m_p_outcome_given_predictors_and_parameters ) ;
		calculate_coefficients( m_outcomes, m_design_probabilities, m_p_outcome_given_predictors_and_parameters, &m_coefficients ) ;
	}
	
	void LogisticRegressionLogLikelihood::calculate_p_outcome_given_predictors_and_parameters(
		Vector const& parameters,
		Matrix const& design_matrices,
		Matrix* result
	) const {
		result->resize( 2, design_matrices.cols() ) ;
		for( int l = 0; l < design_matrices.cols(); ++l ) {
			double const odds = std::exp( parameters.dot( design_matrices.col( l ) ) ) ;
			(*result)( 0, l ) = 1.0 / ( 1.0 + odds ) ;
			(*result)( 1, l ) = odds / ( 1.0 + odds ) ;
		}
	}

	// Calculate the z coefficients, i.e.
	// z_i = ( \psi_i - p_i ) * p_i
	// where \psi_i is the outcome for sample i, and
	//
	// p_i = p( predictors | outcome_i, parameters, prior information )
	//
	// By Bayes theorem p_i is proportional to
	//
	// p( outcome_i | predictors, parameters, prior information ) p( predictors | parameters, prior information )
	//
	// The assumption is that the second term does not depend on the parameters, i.e. that it reduces
	// to p( predictors | prior information ), given by the design probabilities.
	//
	void LogisticRegressionLogLikelihood::calculate_coefficients(
		Vector const& outcomes,
		Matrix const& design_probabilities,
		Matrix const& p_outcome_given_predictors_and_parameters,
		Matrix* result
	) const {
		assert( result != 0 ) ;
		(*result) = design_probabilities  ;
		assert( result->rows() == outcomes.size() ) ;
		
		// Calculate p_i
		for( int i = 0; i < outcomes.size(); ++i ) {
			result->row( i ) = result->row(i).cwiseProduct( p_outcome_given_predictors_and_parameters.row( int( outcomes( i ) )) ) ;
			result->row( i ) /= result->row( i ).sum() ;
		}
		
		// Calculate z_i
		for( int i = 0; i < outcomes.size(); ++i ) {
			for( std::size_t g = 0; g < 3; ++g ) {
				result->operator()( i, g ) *= ( outcomes( i ) - p_outcome_given_predictors_and_parameters( 1, g ) ) ;
			}
		}
	}

	double LogisticRegressionLogLikelihood::get_value_of_function() const {
		return calculate_function(
			m_outcomes,
			m_design_probabilities,
			m_p_outcome_given_predictors_and_parameters
		) ;
	}

	LogisticRegressionLogLikelihood::Vector LogisticRegressionLogLikelihood::get_value_of_first_derivative() const {
		return calculate_first_derivative(
			m_outcomes,
			m_p_outcome_given_predictors_and_parameters,
			m_coefficients,
			m_design_matrices
		) ;
	}

	LogisticRegressionLogLikelihood::Matrix LogisticRegressionLogLikelihood::get_value_of_second_derivative() const {
		return calculate_second_derivative(
			m_outcomes,
			m_p_outcome_given_predictors_and_parameters,
			m_coefficients,
			m_design_matrices
		) ;
	}

	double LogisticRegressionLogLikelihood::calculate_function(
		Vector const& outcomes,
		Matrix const& design_probabilities,
		Matrix const& p_outcome_given_predictors_and_parameters
	) const {
		assert( design_probabilities.rows() == outcomes.size() ) ;
		double result = 0.0 ;
		for( int i = 0; i < outcomes.size(); ++i ) {
			double x = 0.0 ;
			for( int l = 0; l < design_probabilities.size(); ++l ) {
				x += design_probabilities( i, l ) * p_outcome_given_predictors_and_parameters( int( outcomes( i )), l ) ;
			}
			result += std::log( x ) ;
		}
		return result ;
	}

	LogisticRegressionLogLikelihood::Vector LogisticRegressionLogLikelihood::calculate_first_derivative(
		Vector const& outcomes,
		Matrix const& p_outcome_given_predictors_and_parameters,
		Matrix const& coefficients,
		Matrix const& design_matrices
	) const {
		Vector result = Vector::Zero( design_matrices.rows() ) ;
		for( int i = 0; i < outcomes.size(); ++i ) {
			for( int l = 0; l < design_matrices.cols(); ++l ) {
				result += coefficients( i, l ) * design_matrices.col( l ) ;
			}
		}
		return result ;
	}

	LogisticRegressionLogLikelihood::Matrix LogisticRegressionLogLikelihood::calculate_second_derivative(
		Vector const& outcomes,
		Matrix const& p_outcome_given_predictors_and_parameters,
		Matrix const& coefficients,
		Matrix const& design_matrices
	) const {
		assert( coefficients.cols() == design_matrices.cols() ) ;
		Matrix result = Matrix::Zero( design_matrices.rows(), design_matrices.rows() ) ;
		for( int i = 0; i < outcomes.size(); ++i ) {
			Vector second_term = Vector::Zero( design_matrices.rows() ) ;
			for( int l = 0; l < design_matrices.cols(); ++l ) {
				result += coefficients(i, l)
					* ( 1.0 - 2.0 * p_outcome_given_predictors_and_parameters( 1, l ))
					* ( design_matrices.col(l) * design_matrices.col(l).transpose() ) ;
				second_term += coefficients(i, l) * design_matrices.col(l) ;
			}
			result -= second_term * second_term.transpose() ;
		}
		return result ;
	}
	
}
