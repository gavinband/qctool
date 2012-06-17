
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/SNPIdentifyingData.hpp"
#include "components/SNPSummaryComponent/AutosomalFrequentistCaseControlAssociationTest.hpp"

AutosomalFrequentistCaseControlAssociationTest::AutosomalFrequentistCaseControlAssociationTest(
	genfile::CohortIndividualSource const&,
	Vector const& phenotypes,
	Matrix const& covariates
):
	FrequentistCaseControlAssociationTest( phenotypes, covariates )
{}

void AutosomalFrequentistCaseControlAssociationTest::operator()(
	SNPIdentifyingData const& snp,
	Matrix const& genotypes,
	SampleSexes const&, 
	genfile::VariantDataReader&,
	ResultCallback callback
) {
	if( snp.get_position().chromosome().is_autosome() ) {
		return test(
			snp,
			genotypes,
			Vector::LinSpaced( genotypes.cols(), 0, 2 ),
			callback
		) ;
	}
}

