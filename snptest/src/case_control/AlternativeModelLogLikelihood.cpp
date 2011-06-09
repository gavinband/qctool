#include <iostream>
#include <vector>
#include <utility>
#include <boost/iterator/counting_iterator.hpp>
#include "snptest/case_control/AlternativeModelLogLikelihood.hpp"

namespace snptest {
	namespace case_control {
		AlternativeModelLogLikelihood::AlternativeModelLogLikelihood(
			Vector const& phenotypes,
			FinitelySupportedFunctionSet const& genotypes
		):
			LogLikelihood( phenotypes, genotypes )
		{}

		AlternativeModelLogLikelihood::AlternativeModelLogLikelihood(
			Vector const& phenotypes,
			FinitelySupportedFunctionSet const& genotypes,
			std::vector< std::size_t > const& included_samples
		):
			LogLikelihood( phenotypes, genotypes, included_samples )
		{}
		
		double AlternativeModelLogLikelihood::calculate_p_g_theta(
			Vector const& parameters,
			double genotype,
			double phenotype
		) const {
			double x = std::exp( parameters(0) + parameters(1) * genotype ) ;
			return (( phenotype == 1.0 ) ? x : 1.0 ) / ( 1 + x ) ;
		}
	}
}
