#ifndef FREQUENTIST_CASE_CONTROL_ASSOCIATION_TEST_HPP
#define FREQUENTIST_CASE_CONTROL_ASSOCIATION_TEST_HPP

#include <boost/math/distributions/chi_squared.hpp>
#include <Eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "snptest/case_control/LogLikelihood.hpp"
#include "snptest/case_control/NullModelLogLikelihood.hpp"
#include "components/SNPSummaryComponent/AssociationTest.hpp"

struct FrequentistCaseControlAssociationTest: public AssociationTest {
	FrequentistCaseControlAssociationTest(
		Vector const& phenotypes,
		Matrix const& covariates
		) ;

	void operator()( SNPIdentifyingData const& snp, Matrix const& genotypes, genfile::VariantDataReader& data_reader, ResultCallback callback ) ;

private:
	snptest::case_control::LogLikelihood m_alternative_ll ;
	snptest::case_control::NullModelLogLikelihood m_null_ll ;
	Vector const m_phenotypes ;
	Matrix const m_covariates ;
	boost::math::chi_squared_distribution< double > m_chi_squared ;	

private:
	Vector get_genotype_levels( Matrix const& genotypes ) const ;
} ;

#endif
