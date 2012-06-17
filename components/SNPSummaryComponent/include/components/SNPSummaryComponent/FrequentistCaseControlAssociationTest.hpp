
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef FREQUENTIST_CASE_CONTROL_ASSOCIATION_TEST_HPP
#define FREQUENTIST_CASE_CONTROL_ASSOCIATION_TEST_HPP

#include <boost/math/distributions/chi_squared.hpp>
#include <Eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "snptest/case_control/LogLikelihood.hpp"
#include "snptest/case_control/NullModelLogLikelihood.hpp"
#include "components/SNPSummaryComponent/AssociationTest.hpp"

struct FrequentistCaseControlAssociationTest: public AssociationTest {
	FrequentistCaseControlAssociationTest(
		Vector const& phenotypes,
		Matrix const& covariates
	) ;

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const { return prefix + "FrequentistCaseControlAssociationTest" ; }

	void test( SNPIdentifyingData const& snp, Matrix const& predictor_probs, Vector predictor_levels, ResultCallback callback ) ;
	
private:
	snptest::case_control::LogLikelihood m_alternative_ll ;
	snptest::case_control::NullModelLogLikelihood m_null_ll ;
	
	Vector const m_phenotypes ;
	Matrix const m_covariates ;
	boost::math::chi_squared_distribution< double > m_chi_squared ;	

private:
	Vector mean_centre_predictor_levels( genfile::SNPIdentifyingData const& snp, Matrix const& probs, Vector const& predictor_levels ) const ;
} ;

#endif
