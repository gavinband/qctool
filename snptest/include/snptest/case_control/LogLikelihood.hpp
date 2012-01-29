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
			
			LogLikelihood() ;

			LogLikelihood& set_phenotypes( Vector const& phenotypes ) ;
			LogLikelihood& set_genotypes( Matrix const& genotypes, Vector const& levels ) ;
			LogLikelihood& set_covariates( Matrix const& covariates ) ;

			void evaluate_at( Point const& parameters ) ;
			double get_value_of_function() const ;
			Vector get_value_of_first_derivative() const ;
			Matrix get_value_of_second_derivative() const ;

		private:
			Vector m_phenotypes ;
			Matrix m_covariates ;
			Matrix m_genotype_call_probabilities ;
			Vector m_genotype_levels ;
			std::vector< int > m_exclusions ;

			Matrix m_design_matrix ;
			Matrix m_design_matrix_row_tensor_squares ;
			Matrix m_outcome_probabilities ;
			Matrix m_V ;

			double m_value_of_function ;
			Vector m_value_of_first_derivative ;
			Matrix m_value_of_second_derivative ;

		private:
			void add_exclusions( Vector const& v ) ;
			void add_exclusions( Matrix const& v ) ;
			Matrix calculate_design_matrix( Matrix const& covariates ) const ;
			void compute_tensor_squares_of_rows( Matrix& design_matrix, Matrix& result ) ;
			// Calculate the probability of outcome given the genotype, parameters, and covariates.
			Vector evaluate_mean_function( Vector const& linear_combinations, Vector const& outcomes ) const ;
			// Calculate matrix of probabilities of outcome per genotype, given the parameters.
			void calculate_outcome_probabilities( Vector const& parameters, Vector const& phenotypes, Matrix& design_matrix, Matrix* result ) const ;
			void compute_value_of_function( Matrix const& V ) ;
			void compute_value_of_first_derivative( Matrix& V ) ;
			void compute_value_of_second_derivative( Matrix const& V ) ;
		} ;
	}
}

#endif
