#ifndef SNPTEST2_ALTERNATIVE_MODEL_HPP
#define SNPTEST2_ALTERNATIVE_MODEL_HPP

#include <vector>
#include "Eigen/Core"
#include "snptest/case_control/LogLikelihood.hpp"

namespace snptest {
	namespace case_control {
		struct AlternativeModelLogLikelihood: public LogLikelihood
		{
		public:
			typedef Eigen::VectorXd Point ;
			typedef Eigen::VectorXd Vector ;
			typedef Eigen::MatrixXd Matrix ;
			
			AlternativeModelLogLikelihood(
				Vector const& phenotypes,
				Matrix const& genotypes,
				std::vector< double > const& genotype_levels
			) ;

			AlternativeModelLogLikelihood(
				Vector const& phenotypes,
				Matrix const& genotypes,
				std::vector< double > const& genotype_levels,
				std::vector< int > const& included_samples
			) ;

		private:
			double calculate_p_g_theta( Vector const& parameters, double genotype, double phenotype ) const ;
		} ;
	}
}

#endif
