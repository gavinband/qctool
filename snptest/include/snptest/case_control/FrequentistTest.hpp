
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST_CASECONTROL_FREQUENTISTTEST_HPP
#define SNPTEST_CASECONTROL_FREQUENTISTTEST_HPP

#include <boost/math/distributions/chi_squared.hpp>
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "snptest/PerSnpFrequentistTest.hpp"
#include "snptest/FinitelySupportedFunctionSet.hpp"
#include "snptest/case_control/AlternativeModelLogLikelihood.hpp"
#include "snptest/case_control/NullModelLogLikelihood.hpp"

namespace snptest {
	namespace case_control {
		struct FrequentistTest: public PerSnpFrequentistTest {
			FrequentistTest(
				double sample_inclusion_threshhold = 0.1,
				bool mimic_snptest = false
			) ;
			
			Results test(
				genfile::SNPIdentifyingData const& snp,
				Vector const& phenotype_values,
				FinitelySupportedFunctionSet const& all_genotypes,
				Matrix const& covariate_values,
				std::vector< std::size_t > const& indices_of_samples_to_exclude = std::vector< std::size_t >()
			) const ;
		
		private:
			double m_sample_inclusion_threshold ;
			bool m_mimic_snptest ;
			boost::math::chi_squared_distribution< double > m_chi_squared ;	
			
		} ;	
	}
}

#endif
