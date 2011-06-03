#ifndef SNPTEST_LOGISTIC_REGRESSION_MODEL_HPP
#define SNPTEST_LOGISTIC_REGRESSION_MODEL_HPP

#include <vector>
#include "Eigen/Core"
#include "snptest/FinitelySupportedFunctionSet.hpp"

namespace snptest2 {
	// struct LogisticRegressionLogLikelihood
	// This implements a logistic regression model loglikelihood
	// for known (certain) outcome variables and predictor variables
	// that are potentially uncertain.
	struct LogisticRegressionLogLikelihood
	{
	public:
		typedef snptest::FinitelySupportedFunctionSet FinitelySupportedFunctionSet ;
		typedef Eigen::VectorXd Point ;
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::MatrixXd Matrix ;
		
		LogisticRegressionLogLikelihood( LogisticRegressionLogLikelihood const& other ) ;
		LogisticRegressionLogLikelihood(
			Vector const& outcomes,
			Matrix const& design_matrices,
			Matrix const& design_probabilities
		) ;

		void evaluate_at( Point const& parameters ) ;
		double get_value_of_function() const ;
		Vector get_value_of_first_derivative() const ;
		Matrix get_value_of_second_derivative() const ;

	private:
		Vector const& m_outcomes ;
		// Design matrices.
		// Each column is ( 1 l1 l2 l3 l4 )^t where
		// li is a level of the ith predictor variable.
		// There is one column for each possible choice of these levels.
		Matrix const& m_design_matrices ;
		// Design probabilities.
		// Entry in row i and column j is the probability that the levels
		// for sample i are those given in the design matrix for the same column.
		// (So rows should sum to 1.)
		Matrix const& m_design_probabilities ;
		// Probability of outcome given predictors and parameters.
		Matrix m_p_outcome_given_predictors_and_parameters ;
		// Probability of outcome given predictors and parameters
		Matrix m_coefficients ;

	private:
		#if 0
		// Calculate the design matrices (one for each possible level of the predictors)
		// and their probabilities for each sample.
		// The results is a (P+1) x L matrix whose cols are the design matrix
		// where L is the number of predictor levels and P the number of predictor variables,
		// and a N x L matrix of probabilities, one for each sample and each predictor level.
		void calculate_design(
			std::vector< FinitelySupportedFunctionSet > const& predictors,
			Matrix* design_matrices,
			Matrix* design_probabilities
		) const ;
		#endif
		
		// Calculate  p(outcome| predictor, parameters) for each of the two possible outcomes and each predictor level.
		// The result is a 2 x L matrix, stored in the given argument which must be non-null.
		void calculate_p_outcome_given_predictors_and_parameters(
			Vector const& parameters,
			Matrix const& design_matrices,
			Matrix* result
		) const ;

		// Calculate the z coefficients for each sample and each value of the predictors.
		//
		// This is
		//
		// z(i,predictors) = ( \psi_i - p(i,predictors) ) * p(i,predictors)
		//
		// where \psi_i is the outcome for sample i, and
		//
		// p(i,predictors) = p( predictors | outcome_i, parameters, prior information )
		//
		// The result is a N x L matrix (L is the number of predictor levels) and is stored in the given
		// argument which must be non-null.
		//
		// By Bayes theorem p(i,predictors) is proportional to
		//
		// p( predictors | parameters, prior information ) p( outcome_i | predictors, parameters, prior information )
		//
		// The first and second terms are assumed to have been calculated already and supplied in the
		// design_probabilities and p_outcome_given_predictors_and_parameters arguments.
		//
		void calculate_coefficients(
			Vector const& outcomes,
			Matrix const& design_probabilities,
			Matrix const& p_outcome_given_predictors_and_parameters,
			Matrix* result
		) const ;

		double calculate_function(
			Vector const& outcomes,
			Matrix const& design_probabilities,
			Matrix const& p_outcome_given_predictors_and_parameters
		) const ;

		Vector calculate_first_derivative(
			Vector const& phenotypes,
			Matrix const& p_outcome_given_predictors_and_parameters,
			Matrix const& coefficients,
			Matrix const& design_matrices
		) const ;

		Matrix calculate_second_derivative(
			Vector const& outcomes,
			Matrix const& p_outcome_given_predictors_and_parameters,
			Matrix const& coefficientss,
			Matrix const& design_matrices
		) const ;
	
	} ;
}

#endif
