#ifndef SNPTEST_CASE_CONTROL_LOG_LIKELIHOOD_HPP
#define SNPTEST_CASE_CONTROL_LOG_LIKELIHOOD_HPP

#include <vector>
#include <boost/noncopyable.hpp>
#include "Eigen/Core"
#include "snptest/FinitelySupportedFunctionSet.hpp"

namespace snptest {
	namespace case_control {
		struct LogLikelihood: public boost::noncopyable
		{
		public:
			typedef Eigen::VectorXd Point ;
			typedef Eigen::VectorXd Vector ;
			typedef Eigen::RowVectorXd RowVector ;
			typedef Eigen::MatrixXd Matrix ;
			
			LogLikelihood(
				Vector const& phenotypes,
				FinitelySupportedFunctionSet const& genotypes,
				Matrix const& covariates
			) ;

			LogLikelihood(
				Vector const& phenotypes,
				FinitelySupportedFunctionSet const& genotypes,
				Matrix const& covariates,
				std::vector< std::size_t > const& excluded_samples
			) ;

			void evaluate_at( Point const& parameters ) ;
			double get_value_of_function() const ;
			Vector get_value_of_first_derivative() const ;
			Matrix get_value_of_second_derivative() const ;

		private:
			Vector m_phenotypes ;
			Matrix m_genotype_call_probabilities ;
			Matrix m_genotype_levels ;
			Matrix const& m_covariates ;
			std::vector< std::size_t > m_excluded_samples ;

			Matrix m_outcome_probabilities ;
			Matrix m_design_matrix ;
			Matrix m_coefficients ;
			
			double m_value_of_function ;
			Vector m_value_of_first_derivative ;
			Matrix m_value_of_second_derivative ;

		private:
			void deal_with_exclusions( std::vector< std::size_t > exclusions ) ;
			Matrix calculate_design_matrix( Matrix const& covariates ) const ;
			// Calculate the probability of outcome given the genotype, parameters, and covariates.
			Vector evaluate_mean_function( Vector const& linear_combinations, Vector const& outcomes ) const ;
			// Calculate matrix of probabilities of outcome per genotype, given the parameters.
			Matrix calculate_outcome_probabilities( Vector const& parameters, Vector const& phenotypes, Matrix& design_matrix ) const ;
		} ;
	}
}

#endif
