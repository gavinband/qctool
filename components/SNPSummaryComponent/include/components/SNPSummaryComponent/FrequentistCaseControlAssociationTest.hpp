#ifndef FREQUENTIST_CASE_CONTROL_ASSOCIATION_TEST_HPP
#define FREQUENTIST_CASE_CONTROL_ASSOCIATION_TEST_HPP

#include <boost/math/distributions/chi_squared.hpp>
#include <Eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "snptest/case_control/LogLikelihood.hpp"
#include "snptest/case_control/NullModelLogLikelihood.hpp"
#include "integration/NewtonRaphson.hpp"
#include "integration/Derivative.hpp"
#include "integration/maximisation.hpp"
#include "AssociationTest.hpp"

struct FrequentistCaseControlAssociationTest: public AssociationTest {
	FrequentistCaseControlAssociationTest(
		Vector const& phenotypes,
		Matrix const& covariates
	):
		m_chi_squared( 1.0 )
	{
		m_alternative_ll.set_phenotypes( phenotypes ) ;
		m_alternative_ll.set_covariates( covariates ) ;
		m_null_ll.set_phenotypes( phenotypes ) ;
		m_null_ll.set_covariates( covariates ) ;
	}

	void operator()( SNPIdentifyingData const& snp, Matrix const& genotypes, ResultCallback callback ) {
		using integration::derivative ;
		
		Vector genotype_levels = get_genotype_levels( genotypes ) ;
		
		m_null_ll.set_genotypes( genotypes, genotype_levels ) ;
		m_alternative_ll.set_genotypes( genotypes, genotype_levels ) ;

		Vector null_parameters = integration::maximise_by_newton_raphson( m_null_ll, Vector( Vector::Zero( 1 ) ), 0.00001 ) ;
		Vector alternative_parameters = integration::maximise_by_newton_raphson( m_alternative_ll, Vector( Vector::Zero( 2 )), 0.00001 ) ;

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

	Vector get_genotype_levels( Matrix const& genotypes ) const {
		Vector result( Vector::Zero( 3 )) ;
		// compute average genotype
		double const mean_genotype = ( genotypes.col(1) + 2 * genotypes.col(2) ).sum() / ( 2.0 * genotypes.sum() ) ;
		for( int g = 0; g < 3; ++g ) {
			result( g ) = g - mean_genotype ;
		}
		return result ;
	}

private:
	snptest::case_control::LogLikelihood m_alternative_ll ;
	snptest::case_control::NullModelLogLikelihood m_null_ll ;
	boost::math::chi_squared_distribution< double > m_chi_squared ;	
} ;

#endif
