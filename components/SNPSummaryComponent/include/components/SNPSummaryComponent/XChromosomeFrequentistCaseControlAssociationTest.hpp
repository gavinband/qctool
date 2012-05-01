
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

	void operator()( SNPIdentifyingData const& snp, Matrix const& genotypes, SampleSexes const&, genfile::VariantDataReader& data_reader, ResultCallback callback ) ;

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const { return prefix + "XChromosomeFrequentistCaseControlAssociationTest" ; }

	void set_model( std::string const& model ) ;
	
private:
	genfile::CohortIndividualSource const& m_samples ;
	std::vector< char > m_sexes ;
	typedef std::map< char, std::vector< int > > SampleIndices ;
	SampleIndices const m_samples_by_sex ;
	bool const m_with_X_inactivation ;
	bool m_single_organ_model ;
	
private:
	std::vector< char > get_sexes( genfile::CohortIndividualSource const& samples ) const ;
	std::map< char, std::vector< int > > get_samples_by_sex( std::vector< char > const& sexes ) const ;
	
} ;

#endif
