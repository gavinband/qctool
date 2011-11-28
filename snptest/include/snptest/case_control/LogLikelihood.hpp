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
				std::vector< std::size_t > const& included_samples
			) ;

			void evaluate_at( Point const& parameters ) ;
			double get_value_of_function() const ;
			Vector get_value_of_first_derivative() const ;
			Matrix get_value_of_second_derivative() const ;

		private:
			Vector const& m_phenotypes ;
			FinitelySupportedFunctionSet const& m_genotypes ;
			Matrix const& m_covariates ;
			std::vector< Vector > m_design_matrices ;
			std::vector< std::size_t > m_included_samples ;
			Matrix m_p_g_thetas ;
			Matrix m_coefficients ;

		private:
			std::vector< Vector > calculate_design_matrices( Vector const& genotype_levels ) const ;
			virtual double calculate_p_g_theta( Vector const& parameters, double genotype, double phenotype ) const = 0 ;
			Matrix calculate_p_g_thetas( Vector const& parameters, std::vector< Vector > const& design_matrices ) const ;
			Matrix calculate_coefficients( Matrix const& genotypes, Vector const& phenotypes, Matrix const& p_g_thetas ) const ;
			double calculate_function(
				Vector const& phenotypes,
				Matrix const& genotypes,
				Matrix const& p_g_thetas
			) const ;
			Vector calculate_first_derivative(
				Vector const& phenotypes,
				Matrix const& p_g_thetas,
				Matrix const& zs,
				std::vector< Vector > const& design_matrices
			) const ;
			Matrix calculate_second_derivative(
				Vector const& phenotypes,
				Matrix const& p_g_thetas,
				Matrix const& zs,
				std::vector< Vector > const& design_matrices
			) const ;
		} ;
	}
}

#endif
