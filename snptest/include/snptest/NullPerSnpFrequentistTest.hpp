#ifndef SNPTEST_NULL_PER_SNP_FREQUENTIST_TEST_HPP
#define SNPTEST_NULL_PER_SNP_FREQUENTIST_TEST_HPP

#include <vector>
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "snptest/PerSnpFrequentistTest.hpp"

namespace snptest {
	struct NullPerSnpFrequentistTest: public PerSnpFrequentistTest {
		Results test(
			genfile::SNPIdentifyingData const& snp,
			Vector const& phenotype_values,
			FinitelySupportedFunctionSet const& all_genotypes,
			Matrix const& covariate_values,
			std::vector< std::size_t > const& indices_of_samples_to_exclude = std::vector< std::size_t >()
		) const {
			return Results() ;
		}
	} ;
}

#endif

