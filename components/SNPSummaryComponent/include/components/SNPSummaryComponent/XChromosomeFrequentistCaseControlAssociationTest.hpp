
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef XCHROMOSOME_FREQUENTIST_CASE_CONTROL_ASSOCIATION_TEST_HPP
#define XCHROMOSOME_FREQUENTIST_CASE_CONTROL_ASSOCIATION_TEST_HPP

#include <vector>
#include <map>
#include <Eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "components/SNPSummaryComponent/FrequentistCaseControlAssociationTest.hpp"

struct XChromosomeFrequentistCaseControlAssociationTest: public FrequentistCaseControlAssociationTest {
	XChromosomeFrequentistCaseControlAssociationTest(
		genfile::CohortIndividualSource const& samples,
		Vector const& phenotypes,
		Matrix const& covariates,
		bool with_X_inactivation
	) ;

	void operator()( SNPIdentifyingData const& snp, Matrix const& genotypes, genfile::VariantDataReader& data_reader, ResultCallback callback ) ;

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const { return prefix + "XChromosomeFrequentistCaseControlAssociationTest" ; }
private:
	genfile::CohortIndividualSource const& m_samples ;
	std::vector< char > m_sexes ;
	bool const m_with_X_inactivation ;
	typedef std::map< char, std::vector< int > > SampleIndices ;
	SampleIndices m_samples_by_sex ;
	
private:
	std::vector< char > get_sexes( genfile::CohortIndividualSource const& samples ) const ;
	
	int determine_male_coding_column(
		SNPIdentifyingData const& snp,
		Matrix const& genotypes
	) const ;
} ;

#endif
