#include <boost/math/distributions/chi_squared.hpp>
#include <Eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "snptest/case_control/LogLikelihood.hpp"
#include "snptest/case_control/NullModelLogLikelihood.hpp"
#include "integration/NewtonRaphson.hpp"
#include "integration/Derivative.hpp"
#include "integration/maximisation.hpp"
#include "components/SNPSummaryComponent/AssociationTest.hpp"
#include "components/SNPSummaryComponent/FrequentistCaseControlAssociationTest.hpp"

FrequentistCaseControlAssociationTest::FrequentistCaseControlAssociationTest(
	Vector const& phenotypes,
	Matrix const& covariates
):
	m_phenotypes( phenotypes ),
	m_covariates( covariates ),
	m_chi_squared( 1.0 )
{
	m_alternative_ll.set_phenotypes( phenotypes ) ;
	m_alternative_ll.set_covariates( covariates ) ;
	m_null_ll.set_phenotypes( phenotypes ) ;
	m_null_ll.set_covariates( covariates ) ;
}

namespace impl {
	template< typename LL >
	struct AssociationTestStoppingCondition {
		AssociationTestStoppingCondition( LL const& ll, double tolerance ):
			m_ll( ll ),
			m_current_ll_value( -std::numeric_limits< double >::infinity() ),
			m_tolerance( tolerance ),
			m_iteration( 0 )
		{}
		
		bool operator()(
			FrequentistCaseControlAssociationTest::Vector const& value_of_function
		) {
			double const previous_value = m_current_ll_value ;
			m_current_ll_value = m_ll.get_value_of_function() ;
			//std::cerr << "Iteration: " << m_iteration++ << ", previous value: " << previous_value << ", current_value: " << m_current_ll_value << ", magnitude of difference: " << std::abs( m_current_ll_value - previous_value ) << ".\n" ;
			return diverged() || ( std::abs( m_current_ll_value - previous_value ) < m_tolerance ) ;
		}

		bool diverged() const { return m_current_ll_value != m_current_ll_value ; }
		std::size_t number_of_iterations() const { return m_iteration ; }

		private:
			LL const& m_ll ;
			double m_current_ll_value ;
			double const m_tolerance ;
			std::size_t m_iteration ;
			bool m_converged ;
	} ;
}

void FrequentistCaseControlAssociationTest::operator()( SNPIdentifyingData const& snp, Matrix const& genotypes, genfile::VariantDataReader&, ResultCallback callback ) {
	using integration::derivative ;
	
	Vector genotype_levels = get_genotype_levels( genotypes ) ;
	
	m_null_ll.set_genotypes( genotypes, genotype_levels ) ;
	m_alternative_ll.set_genotypes( genotypes, genotype_levels ) ;

	Vector null_parameters = Vector::Zero( m_covariates.cols() + 1 ) ;
	{
		impl::AssociationTestStoppingCondition< snptest::case_control::NullModelLogLikelihood > stopping_condition( m_null_ll, 0.01 ) ;
		null_parameters = integration::maximise_by_newton_raphson(
			m_null_ll,
			null_parameters,
			stopping_condition
		) ;
	}
	Vector alternative_parameters = Vector::Zero( null_parameters.size() + 1 ) ;
	alternative_parameters( 0 ) = null_parameters(0) ;
	alternative_parameters.tail( m_covariates.cols() ) = null_parameters.tail( m_covariates.cols() ) ;
	{
		impl::AssociationTestStoppingCondition< snptest::case_control::LogLikelihood > stopping_condition( m_alternative_ll, 0.01 ) ;
		alternative_parameters = integration::maximise_by_newton_raphson(
			m_alternative_ll,
			alternative_parameters,
			stopping_condition
		) ;
	}
	callback( "null_ll", m_null_ll.get_value_of_function() ) ;
	callback( "alternative_ll", m_alternative_ll.get_value_of_function() ) ;
	double const test_statistic = -2.0 * ( m_null_ll.get_value_of_function() - m_alternative_ll.get_value_of_function() ) ;
	callback( "test_statistic", test_statistic ) ;
	if( test_statistic > 0.0 ) {
		double const p_value = boost::math::cdf(
			boost::math::complement(
				m_chi_squared,
				test_statistic
			)
		) ;
		callback( "p_value", p_value ) ;
	}
	
	using genfile::string_utils::to_string ;
	for( int i = 0; i < alternative_parameters.size(); ++i ) {
		callback( "beta_" + to_string( i ), alternative_parameters( i ) ) ;
	}

	Matrix variance_covariance = ( -m_alternative_ll.get_value_of_second_derivative() ).inverse() ;
	for( int i = 0; i < variance_covariance.rows(); ++i ) {
		for( int j = 0; j < variance_covariance.cols(); ++j ) {
			callback( "sigma_" + to_string( i ) + "," + to_string( j ), variance_covariance( i, j ) ) ;
		}
	}
	for( int i = 0; i < variance_covariance.rows(); ++i ) {
		callback( "se_" + to_string( i ), std::sqrt( variance_covariance( i, i )) ) ;
	}
}

FrequentistCaseControlAssociationTest::Vector FrequentistCaseControlAssociationTest::get_genotype_levels( Matrix const& genotypes ) const {
	Vector result( Vector::Zero( 3 )) ;
	// compute average genotype
	double const mean_genotype = ( genotypes.col(1) + 2.0 * genotypes.col(2) ).sum() / genotypes.sum() ;
	result <<
		-mean_genotype, 1.0 - mean_genotype, 2.0 - mean_genotype ;
	return result ;
}
