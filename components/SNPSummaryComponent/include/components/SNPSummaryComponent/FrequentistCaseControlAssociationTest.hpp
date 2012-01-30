#ifndef FREQUENTIST_CASE_CONTROL_ASSOCIATION_TEST_HPP
#define FREQUENTIST_CASE_CONTROL_ASSOCIATION_TEST_HPP

#include "AssociationTest.hpp"
#include "snptest/case_control/LogLikelihood.hpp"
#include "snptest/case_control/NullModelLogLikelihood.hpp"
#include "integration/NewtonRaphson.hpp"

struct FrequentistCaseControlAssociationTest: public AssociationTest {
	FrequentistCaseControlAssociationTest(
		Vector const& phenotypes,
		Matrix const& covariates
	) {
		m_alternative_ll.set_phenotypes( phenotypes ) ;
		m_alternative_ll.set_covariates( covariates ) ;
		m_null_ll.set_phenotypes( phenotypes ) ;
		m_null_ll.set_covariates( covariates ) ;
	}

	void operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, ResultCallback ) const {
		using integration::derivative ;
		
		m_null_loglikelihood.set_genotypes( genotypes ) ;
		integration::Derivative< snptest::case_control::NullModelLogLikelihood >
			null_loglikelihood_derivative = derivative( m_null_loglikelihood ) ;
		Vector null_parameters = integration::find_root_by_newton_raphson(
			null_loglikelihood_derivative,
			Vector::Zero( 1 ),
			0.00001
		) ;

		m_alternative_loglikelihood.set_genotypes( genotypes ) ;
		integration::Derivative< snptest::case_control::AlternativeModelLogLikelihood >
			alternative_loglikelihood_derivative = derivative( m_alternative_loglikelihood ) ;
		Vector alternative_parameters = integration::find_root_by_newton_raphson(
			alternative_loglikelihood_derivative,
			Vector::Zero( 2 ),
			0.00001
		) ;

		callback( "null_loglikelihood", m_null_loglikelihood.get_value_of_function() ) ;
		callback( "alternative_loglikelihood", m_alternative_loglikelihood.get_value_of_function() ) ;
		double const test_statistic = -2.0 * ( m_null_loglikelihood.get_value_of_function() - m_alternative_loglikelihood.get_value_of_function() ;
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

		Matrix variance_covariance = (-alternative_loglikelihood.get_value_of_second_derivative()).inverse() ;
		for( int i = 0; i < variance_covariance.rows(); ++i ) {
			for( int j = 0; j < variance_covariance.cols(); ++j ) {
				callback( "sigma_" + to_string( i ) + "," + to_string( j ), variance_covariance( i, j ) ) ;
			}
		}
		for( int i = 0; i < variance_covariance.rows(); ++i ) {
			callback( "se_" + to_string( i ), std::sqrt( variance_covariance( i, i )) ) ;
		}
	}

private:
	snptest::case_control::LogLikelihood m_alternative_ll ;
	snptest::case_control::NullModelLogLikelihood m_null_ll ;
} ;

#endif
