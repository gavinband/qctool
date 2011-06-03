#include <boost/math/distributions.hpp>
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/Error.hpp"
#include "integration/NewtonRaphson.hpp"
#include "integration/Derivative.hpp"
#include "snptest/case_control/FrequentistTest.hpp"
#include "snptest/case_control/NullModelLogLikelihood.hpp"
#include "snptest/case_control/AlternativeModelLogLikelihood.hpp"

namespace snptest {
	namespace case_control {
		FrequentistTest::FrequentistTest():
			m_fill_null_genotypes( false ),
			m_chi_squared( 1.0 )
		{}
	
		FrequentistTest::Results FrequentistTest::test(
			Vector const& phenotypes,
			Matrix const& covariates,	
			genfile::SingleSNPGenotypeProbabilities const& all_genotypes,
			std::vector< std::size_t > const& indices_of_samples_to_include
		) const {
			assert( std::size_t( phenotypes.size() ) == indices_of_samples_to_include.size() ) ;
			if( covariates.rows() != 0 || covariates.cols() != 0 ) {
				throw genfile::BadArgumentError( "CaseControlFrequentistTest::test()", "covariates (nonempty)" ) ;
			}

			Matrix genotypes( phenotypes.size(), 3 ) ;
			std::vector< double > genotype_levels(3) ;
	
			{
				for( int i = 0; i < phenotypes.size(); ++i ) {
					for( std::size_t g = 0; g < 3; ++g ) {
						genotypes( i, g ) = all_genotypes( indices_of_samples_to_include[i], g ) ;
					}
				}
				// We now deal with missing genotypes (= null calls).
				// We fill in any null calls from the estimated allele frequency.
				// The allele frequency estimate is a maximum likelihood estimate,
				// but for the likelihood which treates the calls as the actual data.
				// In practice I guess this is ok.
				double I = genotypes.sum() ;
				double theta = ( genotypes.col(1) + 2.0 * genotypes.col(2) ).sum() / (2.0 * I ) ;
				std::size_t const N = genotypes.rows() ;

				Vector expected_genotypes( 3 ) ;
				expected_genotypes
					<< ( 1.0 - theta ) * ( 1.0 - theta ),
					2.0 * theta * ( 1.0 - theta ),
					theta * theta ;
				if( m_fill_null_genotypes ) {
					for( std::size_t i = 0; i < N; ++i ) {
						double sum = genotypes.row(i).sum() ;
						if( sum > 1.0 ) {
							// assume this is due to rounding error.
							genotypes.row(i) /= sum ;
						}
						else {
							// fill in with expected genotype.
							genotypes.row(i) += ( 1.0 - sum ) * expected_genotypes ;
						}
					}
				}
				double mean_genotype = ( genotypes.col(1) + 2.0 * genotypes.col(2) ).sum() / N ;
				for( std::size_t g = 0; g < 3; ++g ) {
					genotype_levels[g] = double( g ) - mean_genotype ;
				}
			}
	
			using integration::derivative ;
			snptest::case_control::NullModelLogLikelihood null_loglikelihood( phenotypes ) ;
			integration::Derivative< snptest::case_control::NullModelLogLikelihood > null_loglikelihood_derivative= derivative( null_loglikelihood ) ;
			Vector null_parameters = integration::find_root_by_newton_raphson(
				null_loglikelihood_derivative,
				Vector::Zero( 1 ),
				0.00001
			) ;

			/*
			m_ui_context.logger() << "AssociationTester:        null: MLE is " << null_parameters << ".\n" ;
			m_ui_context.logger() << "AssociationTester:        null: loglikelihood is " << null_loglikelihood.get_value_of_function() << ".\n" ;

			m_ui_context.logger() << "AssociationTester: alternative: Finding MLE...\n" ;
			*/
			snptest::case_control::AlternativeModelLogLikelihood alternative_loglikelihood(
				phenotypes,
				genotypes,
				genotype_levels
			) ;
			integration::Derivative< snptest::case_control::AlternativeModelLogLikelihood > alternative_loglikelihood_derivative = derivative( alternative_loglikelihood ) ;
			Vector alternative_parameters = integration::find_root_by_newton_raphson(
				alternative_loglikelihood_derivative,
				Vector::Zero( 2 ),
				0.00001
			) ;
/*
			m_ui_context.logger() << "AssociationTester: alternative: MLE is " << alternative_parameters << ".\n" ;
			m_ui_context.logger() << "AssociationTester: alternative: loglikelihood is " << alternative_loglikelihood.get_value_of_function() << ".\n" ;
			m_ui_context.logger() << "AssociationTester: -2LR = " << -2.0 * ( null_loglikelihood.get_value_of_function() - alternative_loglikelihood.get_value_of_function() ) << ".\n" ;
*/
			Results result ;
			result.test_statistic = result.p_value = result.beta = result.standard_error = std::numeric_limits< double >::quiet_NaN() ;
			result.test_statistic = -2.0 * ( null_loglikelihood.get_value_of_function() - alternative_loglikelihood.get_value_of_function() ) ;
			if( result.test_statistic > 0.0 ) {
				result.p_value = boost::math::cdf(
					boost::math::complement(
						m_chi_squared,
						-2.0 * ( null_loglikelihood.get_value_of_function() - alternative_loglikelihood.get_value_of_function() )
					)
				) ;
			}
			result.beta = alternative_parameters( 1 ) ;
			result.variance_covariance = (-alternative_loglikelihood.get_value_of_second_derivative()).inverse() ;
			result.standard_error = std::sqrt( result.variance_covariance( 1, 1 ) ) ;
			return result ;
		}
	}
}

