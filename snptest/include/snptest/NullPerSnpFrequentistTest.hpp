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
			genfile::SingleSNPGenotypeProbabilities const& all_genotypes,
			std::vector< std::size_t > const& indices_of_samples_to_include
		) const {
			return Results() ;
		}
	} ;
}

#endif

