
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST2_NULL_MODEL_HPP
#define SNPTEST2_NULL_MODEL_HPP

#include <boost/noncopyable.hpp>
#include <vector>
#include "Eigen/Eigen"
#include "snptest/FinitelySupportedFunctionSet.hpp"

namespace snptest {
	namespace case_control {
		struct NullModelLogLikelihood: public boost::noncopyable
		{
		public:
			typedef Eigen::VectorXd Vector ;
			typedef Eigen::RowVectorXd RowVector ;
			typedef Eigen::MatrixXd Matrix ;
			
			NullModelLogLikelihood() ;

			NullModelLogLikelihood& set_phenotypes( Vector const& phenotypes ) ;
			NullModelLogLikelihood& set_covariates( Matrix const& covariates ) ;
			NullModelLogLikelihood& set_predictor_probs( Matrix const& genotypes, Vector const& levels ) ;
			NullModelLogLikelihood& add_exclusions( std::vector< int > const& exclusions ) ;

			void evaluate_at( Vector const& parameters ) ;
			double get_value_of_function() const ;
			Vector get_value_of_first_derivative() const ;
			Matrix get_value_of_second_derivative() const ;

		private:
			Vector m_phenotypes ;
			Matrix m_covariates ;
			Matrix m_genotype_call_probabilities ;
			Vector m_genotype_levels ;
			std::vector< int > m_exclusions ;

			Matrix m_outcome_probabilities ;
			Matrix m_design_matrix ;
			Matrix m_coefficients ;
			
			double m_value_of_function ;
			Vector m_value_of_first_derivative ;
			Matrix m_value_of_second_derivative ;

		private:
			void add_exclusions( Vector const& v ) ;
			void add_exclusions( Matrix const& matrix ) ;
			
			Matrix calculate_design_matrix( Matrix const& covariates ) const ;
			// Calculate the probability of outcome given the genotype, parameters, and covariates.
			Vector evaluate_mean_function( Vector const& linear_combinations, Vector const& outcomes ) const ;
			// Calculate matrix of probabilities of outcome per genotype, given the parameters.
			Matrix calculate_outcome_probabilities( Vector const& parameters, Vector const& phenotypes, Matrix& design_matrix ) const ;
		} ;
	}
}

#endif
