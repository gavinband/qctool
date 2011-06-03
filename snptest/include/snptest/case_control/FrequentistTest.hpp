#ifndef SNPTEST_CASECONTROL_FREQUENTISTTEST_HPP
#define SNPTEST_CASECONTROL_FREQUENTISTTEST_HPP

#include <boost/math/distributions/chi_squared.hpp>
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "snptest/PerSnpFrequentistTest.hpp"

namespace snptest {
	namespace case_control {
		struct FrequentistTest: public PerSnpFrequentistTest {
			FrequentistTest() ;
		
			Results test(
				Vector const& phenotype_values,
				Matrix const& covariate_values,
				genfile::SingleSNPGenotypeProbabilities const& all_genotypes,
				std::vector< std::size_t > const& indices_of_samples_to_include
			) const ;
		
		private:
			bool m_fill_null_genotypes ;
			boost::math::chi_squared_distribution< double > m_chi_squared ;	
		} ;	
	}
}

#endif
