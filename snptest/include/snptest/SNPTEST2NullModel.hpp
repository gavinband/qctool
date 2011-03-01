#ifndef SNPTEST2_NULL_MODEL_HPP
#define SNPTEST2_NULL_MODEL_HPP

#include <vector>
#include "Eigen/Eigen"

namespace snptest2 {
	
	struct NullModelLogLikelihood
	{
	public:
		typedef Eigen::VectorXd Point ;
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::MatrixXd Matrix ;
			
		NullModelLogLikelihood(
			Vector const& phenotypes
		) ;
		
		void evaluate_at( Vector const& parameters ) ;
		
		double get_value_of_function() const ;
		Vector get_value_of_first_derivative() const ;
		Matrix get_value_of_second_derivative() const ;
		
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
		
		Vector const& m_phenotypes ;
		Matrix m_p_thetas ;
	} ;
	
	typedef NullModelLogLikelihood NullModelEvaluator ;
}

#endif
