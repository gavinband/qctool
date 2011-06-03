#ifndef SNPTEST2_ALTERNATIVE_MODEL_HPP
#define SNPTEST2_ALTERNATIVE_MODEL_HPP

#include <vector>
#include "Eigen/Core"

namespace snptest {
	namespace case_control {
		struct AlternativeModelLogLikelihood
		{
		public:
			typedef Eigen::VectorXd Point ;
			typedef Eigen::VectorXd Vector ;
			typedef Eigen::MatrixXd Matrix ;
			
			AlternativeModelLogLikelihood( AlternativeModelLogLikelihood const& other ) ;
			AlternativeModelLogLikelihood(
				Vector const& phenotypes,
				Matrix const& genotypes,
				std::vector< double > const& genotype_levels
			) ;

			void evaluate_at( Point const& parameters ) ;
			double get_value_of_function() const ;
			Vector get_value_of_first_derivative() const ;
			Matrix get_value_of_second_derivative() const ;

		private:
			std::vector< Vector > calculate_design_matrices( std::vector< double > const& genotype_levels ) const ;
			double calculate_p_g_theta( Vector const& parameters, double genotype, double phenotype ) const ;
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
		
		private:
			Vector const& m_phenotypes ;
			Matrix const& m_genotypes ;
			std::vector< Vector > m_design_matrices ;
			Matrix m_p_g_thetas ;
			Matrix m_coefficients ;
		} ;
	}
}

#endif
