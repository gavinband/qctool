
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef AUTOSOMAL_FREQUENTIST_CASE_CONTROL_ASSOCIATION_TEST_HPP
#define AUTOSOMAL_FREQUENTIST_CASE_CONTROL_ASSOCIATION_TEST_HPP

#include <Eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "components/SNPSummaryComponent/FrequentistCaseControlAssociationTest.hpp"

struct AutosomalFrequentistCaseControlAssociationTest: public FrequentistCaseControlAssociationTest {
	AutosomalFrequentistCaseControlAssociationTest(
		genfile::CohortIndividualSource const&,
		Vector const& phenotypes,
		Matrix const& covariates
	) ;

	void operator()( SNPIdentifyingData const& snp, Matrix const& genotypes, genfile::VariantDataReader& data_reader, ResultCallback callback ) ;

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const { return prefix + "AutosomalFrequentistCaseControlAssociationTest" ; }
} ;

#endif
