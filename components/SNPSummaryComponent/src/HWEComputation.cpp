
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include "genfile/Error.hpp"
#include "metro/likelihood/Multinomial.hpp"
#include "components/SNPSummaryComponent/SNPHWE.hpp"
#include "components/SNPSummaryComponent/HWEComputation.hpp"
#include "metro/likelihood/Multinomial.hpp"
#include "metro/likelihood/ProductOfMultinomials.hpp"
#include "metro/FishersExactTest.hpp"

// #define DEBUG_HWE_COMPUTATION 1

namespace snp_summary_component {
	HWEComputation::HWEComputation():
		m_threshhold( 0.9 ),
		m_chi_squared_1df( 1.0 ),
		m_chi_squared_2df( 2.0 )
	{}
	
	void HWEComputation::operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, SampleSexes const& sexes, genfile::VariantDataReader&, ResultCallback callback ) {
		genfile::Chromosome const& chromosome = snp.get_position().chromosome() ;
		if( chromosome == genfile::Chromosome( "0X" ) ) {
			X_chromosome_test( snp, genotypes, sexes, callback ) ;
		}
		else if( chromosome.is_autosome() ) {
			autosomal_test( snp, genotypes, callback ) ;
		}
	}

	std::string HWEComputation::get_summary( std::string const& prefix, std::size_t column_width ) const { return prefix + "HWEComputation" ; }

	void HWEComputation::autosomal_test( SNPIdentifyingData const& snp, Genotypes const& genotypes, ResultCallback callback ) {
		Eigen::VectorXd genotype_counts = Eigen::VectorXd::Zero( 3 ) ;
		for( int g = 0; g < 3; ++g ) {
			genotype_counts( g ) = std::floor( genotypes.col(g).sum() + 0.5 ) ;
		}
		autosomal_exact_test( snp, genotype_counts, callback ) ;
		autosomal_multinomial_test( snp, genotype_counts, callback ) ;
	}

	void HWEComputation::autosomal_exact_test( SNPIdentifyingData const& snp, Eigen::VectorXd const& genotype_counts, ResultCallback callback ) {
		if( genotype_counts.array().maxCoeff() > 0.5 ) {
			double HWE_pvalue = SNPHWE( genotype_counts(1), genotype_counts(0), genotype_counts(2) ) ;
			callback( "HW_exact_p_value", HWE_pvalue ) ;
		}
		else {
			callback( "HW_exact_p_value", genfile::MissingValue() ) ;
		}
	}

	void HWEComputation::autosomal_multinomial_test( SNPIdentifyingData const& snp, Eigen::VectorXd const& genotype_counts, ResultCallback callback ) {
		metro::likelihood::Multinomial< double, Eigen::VectorXd, Eigen::MatrixXd > hw_model( genotype_counts ) ;
		{
			// compute MLE under assumption of hardy-weinberg.
			double p = 2.0 * genotype_counts( 2 ) + genotype_counts( 1 ) ;
			p /= 2.0 * genotype_counts.sum() ;
			Eigen::VectorXd a = Eigen::VectorXd::Zero( 3 ) ;
			a( 0 ) = ( 1 - p ) * ( 1 - p ) ;
			a( 1 ) = 2.0 * p * ( 1 - p ) ;
			a( 2 ) = p * p ;
			hw_model.evaluate_at( a ) ;
		}

		metro::likelihood::Multinomial< double, Eigen::VectorXd, Eigen::MatrixXd > full_model( genotype_counts ) ;
		full_model.evaluate_at( full_model.get_MLE() ) ;
		
#if DEBUG_HWE_COMPUTATION
		std::cerr << std::fixed << std::setprecision( 5 );
		std::cerr << "table = " << genotype_counts.transpose() << "...\n" ;
		std::cerr << std::resetiosflags( std::ios::floatfield ) ;
		std::cerr << "hw_model.MLE = " << hw_model.get_parameters().transpose() << ", loglikelihood = " << hw_model.get_value_of_function() << "...\n" ;
		std::cerr << "full_model.MLE = " << full_model.get_parameters().transpose() << ", loglikelihood = " << full_model.get_value_of_function() << "...\n" ;
#endif
		
		double likelihood_ratio_statistic = 2.0 * ( full_model.get_value_of_function() - hw_model.get_value_of_function() ) ;
		double p_value = std::numeric_limits< double >::quiet_NaN() ;
		if( likelihood_ratio_statistic != likelihood_ratio_statistic || likelihood_ratio_statistic < 0.0 ) {
			likelihood_ratio_statistic = std::numeric_limits< double >::quiet_NaN() ;
		}
		else {
			p_value = boost::math::cdf(
				boost::math::complement(
					m_chi_squared_1df,
					likelihood_ratio_statistic
				)
			) ;
		}
		
		callback( "HW_multinomial_p_value", p_value ) ;
	}

	void HWEComputation::X_chromosome_test( SNPIdentifyingData const& snp, Genotypes const& genotypes, SampleSexes const& sexes, ResultCallback callback ) {
		// We look at three models and perform two LR tests.
		// model1: full model, males and females may have different frequencies and no assumption of HW in females.  (3 parameters)
		// model2: HWE holds in females, but males and females may have different frequencies. (2 parameters)
		// model3: HWE holds in females and males and females have the same frequency. (1 parameter)
		//
		// test1: test that males and females have the same allele frequency
		// test2: test HWE in females.

		Eigen::MatrixXd genotype_counts = Eigen::MatrixXd::Zero( 2, 3 ) ; // first row is males, second row is females.  Last column will be zero for males.
		Eigen::MatrixXd allele_counts = Eigen::MatrixXd::Zero( 2, 2 ) ; // first row is males, second row is females.

		typedef Eigen::VectorXd Vector ;

		int const MALES = 0 ;
		int const FEMALES = 1 ;

		for( int i = 0; i < genotypes.rows(); ++i ) {
			if( sexes[ i ] == 'm' || sexes[ i ] == 'f' ) {
				int const index = ( sexes[ i ] == 'm' ) ? MALES : FEMALES ;
				for( int g = 0; g < 3; ++g ) {
					if( genotypes( i, g ) > m_threshhold ) {
						++genotype_counts( index, g ) ;
						break ;
					}
				}
			}
		}
		for( int g = 0; g < 2; ++g ) {
			allele_counts( 0, g ) = genotype_counts( 0, g ) ;
			allele_counts( 1, g ) = genotype_counts( 1, 1 ) + 2.0 * genotype_counts( 1, 2.0 * g ) ;
		}

		if( genotype_counts.maxCoeff() > 0 ) {
			typedef metro::likelihood::Multinomial< double, Eigen::VectorXd, Eigen::MatrixXd > Multinomial ;
			typedef metro::likelihood::ProductOfMultinomials< double, Eigen::VectorXd, Eigen::MatrixXd > ProductOfIndependentMultinomials ;
			// In model 1 each genotype is allowed to have its own frequency.
			// Morever males and females may have different frequencies
			ProductOfIndependentMultinomials full_model( genotype_counts ) ;
			full_model.evaluate_at( full_model.get_MLE() ) ;

			// In model2 alleles are independent (HWE) in females, but males and females may differ.
			ProductOfIndependentMultinomials model2( genotype_counts ) ;
			{
				double const p_females = ( 2 * genotype_counts( FEMALES, 0 ) + genotype_counts( FEMALES, 1 ) ) / (2 * genotype_counts.row( FEMALES ).sum() ) ;
				double const p_males = genotype_counts( MALES, 0 ) / genotype_counts.row( MALES ).sum() ;
				Vector params( 6 ) ;
				params( 0 ) = p_males ;
				params( 1 ) = 1-p_males ;
				params( 2 ) = 0 ;
				params( 3 ) = p_females * p_females ;
				params( 4 ) = 2 * p_females * (1-p_females) ;
				params( 5 ) = (1-p_females) * (1-p_females) ;
				model2.evaluate_at( params ) ;
			}

			// In model3, the NULL model, alleles are independent (HWE) and males and females must agree.
			ProductOfIndependentMultinomials model3( genotype_counts ) ;
			{
				double const p = ( 2 * genotype_counts( FEMALES, 0 ) + genotype_counts( MALES, 0 ) + genotype_counts( FEMALES, 1 ) )
					/ ( 2 * genotype_counts.row( FEMALES ).sum() + genotype_counts.row( MALES ).sum() ) ;
				
				Vector params( 6 ) ;
				params( 0 ) = p ;
				params( 1 ) = (1-p) ;
				params( 2 ) = 0 ;
				params( 3 ) = p * p ;
				params( 4 ) = 2 * p * (1-p) ;
				params( 5 ) = (1-p) * (1-p) ;
				model3.evaluate_at( params ) ;
			}

#if DEBUG_HWE_COMPUTATION
			std::cerr << "Genotype counts:\n" << genotype_counts << "\n" ;
			std::cerr << "Allele counts:\n" << allele_counts << "\n" ;
			std::cerr << "Model 1 params: " << full_model.get_parameters().transpose() << ".\n" ;
			std::cerr << "Model 2 params: " << model2.get_parameters().transpose() << ".\n" ;
			std::cerr << "Model 3 params: " << model3.get_parameters().transpose() << ".\n" ;
#endif
			double const NaN = std::numeric_limits< double >::quiet_NaN() ;

			using boost::math::cdf ;
			using boost::math::complement ;

			if( genotype_counts.row( FEMALES ).array().maxCoeff() > 0.5 ) {
				double exact_HWE_pvalue = SNPHWE( genotype_counts( FEMALES, 1 ), genotype_counts( FEMALES, 0 ), genotype_counts( FEMALES, 2 ) ) ;
				callback( "HW_females_exact_pvalue", exact_HWE_pvalue ) ;
			} else {
				callback( "HW_females_exact_pvalue", genfile::MissingValue() ) ;
			}
			
			double const lr_stat_12 = 2.0 * ( full_model.get_value_of_function() - model2.get_value_of_function() ) ;
			double p_value_12 = NaN ;
			if( lr_stat_12 == lr_stat_12 && lr_stat_12 > 0 && lr_stat_12 != std::numeric_limits< double >::infinity() ) {
				p_value_12 = cdf( complement( m_chi_squared_1df, lr_stat_12 ) ) ;
				callback( "HWE_females_lr_pvalue", p_value_12 ) ;
			}

			// Also get exact male/female p-value
			{
				Eigen::Matrix2d A = allele_counts ;
				A(0,0) = std::floor( A(0,0) + 0.5 ) ;
				A(0,1) = std::floor( A(0,1) + 0.5 ) ;
				A(1,0) = std::floor( A(1,0) + 0.5 ) ;
				A(1,1) = std::floor( A(1,1) + 0.5 ) ;
				double const male_female_pvalue = metro::FishersExactTest( A ).get_pvalue( metro::FishersExactTest::eTwoSided ) ;
				if( male_female_pvalue == male_female_pvalue ) {
					callback( "male_female_exact_pvalue", male_female_pvalue ) ;
				} else {
					callback( "male_female_exact_pvalue", genfile::MissingValue() ) ;
				}
			}

			double const lr_stat_23 = 2.0 * ( model2.get_value_of_function() - model3.get_value_of_function() ) ;
			double p_value_23 = NaN ;
			if( lr_stat_23 == lr_stat_23 && lr_stat_23 > 0 && lr_stat_23 != std::numeric_limits< double >::infinity() ) {
				p_value_23 = cdf( complement( m_chi_squared_1df, lr_stat_23 ) ) ;
				callback( "male_female_lr_pvalue", p_value_23 ) ;
			}

			double const lr_stat_13 = 2.0 * ( full_model.get_value_of_function() - model3.get_value_of_function() ) ;
			double p_value_13 = NaN ;
			if( lr_stat_13 == lr_stat_13 && lr_stat_13 > 0 && lr_stat_13 != std::numeric_limits< double >::infinity() ) {
				p_value_13 = cdf( complement( m_chi_squared_2df, lr_stat_13 ) ) ;
				callback( "male_female_and_HWE_2df_lr_pvalue", p_value_13 ) ;
			}
		}
	}
}

