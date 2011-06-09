#ifndef SNPTEST_NULL_PER_SNP_FREQUENTIST_TEST_HPP
#define SNPTEST_NULL_PER_SNP_FREQUENTIST_TEST_HPP

#include <vector>
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "snptest/PerSnpFrequentistTest.hpp"

namespace snptest {
	struct NullPerSnpFrequentistTest: public PerSnpFrequentistTest {
		Results test(
			Vector const& phenotype_values,
			Matrix const& covariate_values,
			genfile::SNPIdentifyingData const& snp,
			FinitelySupportedFunctionSet const& all_genotypes
		) const {
			return Results() ;
		}
	} ;
}

#endif

