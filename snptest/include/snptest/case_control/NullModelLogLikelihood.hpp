#ifndef SNPTEST2_NULL_MODEL_HPP
#define SNPTEST2_NULL_MODEL_HPP

#include <boost/noncopyable.hpp>
#include <vector>
#include "Eigen/Eigen"
#include "snptest/FinitelySupportedFunctionSet.hpp"

namespace snptest {
	namespace case_control {
		struct NullModelLogLikelihood:public boost::noncopyable
		{
		public:
			typedef Eigen::VectorXd Point ;
			typedef Eigen::VectorXd Vector ;
			typedef Eigen::MatrixXd Matrix ;
			
			NullModelLogLikelihood(
				Vector const& phenotypes,
				FinitelySupportedFunctionSet const& genotypes,
				bool weight_by_genotypes = false
			) ;

			NullModelLogLikelihood(
				Vector const& phenotypes,
				FinitelySupportedFunctionSet const& genotypes,
				bool weight_by_genotypes,
				std::vector< std::size_t > const& included_samples
			) ;
		
			void evaluate_at( Vector const& parameters ) ;
		
			double get_value_of_function() const ;
			Vector get_value_of_first_derivative() const ;
			Matrix get_value_of_second_derivative() const ;
		private:
			Vector const& m_phenotypes ;
			FinitelySupportedFunctionSet const& m_genotypes ;
			bool m_weight_by_genotypes ;
			std::vector< std::size_t > const m_included_samples ;
			Matrix m_p_thetas ;
		
		private:
			NullModelLogLikelihood( NullModelLogLikelihood const& other ) ;
		
			std::vector< Vector > calculate_design_matrices() const ;
			double calculate_p_theta( Vector const& parameters, double phenotype ) const ;
			Matrix calculate_p_thetas( Vector const& parameters ) const ;
			double calculate_function(
				Vector const& phenotypes,
				Matrix const& p_thetas
			) const ;
			Vector calculate_first_derivative(
				Vector const& phenotypes,
				Matrix const& p_thetas
			) const ;
			Matrix calculate_second_derivative(
				Vector const& phenotypes,
				Matrix const& p_thetas
			) const ;
		} ;
	
		typedef NullModelLogLikelihood NullModelEvaluator ;
	}
}

#endif
