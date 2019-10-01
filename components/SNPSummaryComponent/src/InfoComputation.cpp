
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <memory>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <Eigen/Core>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/VariantDataReader.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "components/SNPSummaryComponent/InfoComputation.hpp"

namespace stats {
	void InfoComputation::operator()(
		VariantIdentifyingData const& snp,
		Genotypes const& genotypes,
		Ploidy const& ploidy,
		genfile::VariantDataReader&,
		ResultCallback callback
	) {
		m_computation.compute( snp, genotypes, ploidy ) ;
		callback( "info", m_computation.info() ) ;
		callback( "impute_info", m_computation.impute_info() ) ;
	}
	
	std::string InfoComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "InfoComputation" ;
	}
	
	namespace impl {
		namespace {
			double compute_variance( Eigen::VectorXd const& levels, Eigen::VectorXd const& levels_squared, Eigen::VectorXd const& probs ) {
				double const mean = ( probs.transpose() * levels )(0) ;
				double const variance = ( probs.transpose() * levels_squared )(0) - ( mean * mean ) ;
				return variance ;
			}

			// treat the rows of the probabilities matrix as probabilities.
			// distribution for individual i is taken as a mixture of the distribution given by row i of probabilities,
			// and the fallback distribution if the sum of row i < 1, appropriately weighted.
			// Only individuals with inclusion = 1 are used.
			template< typename Probabilities, typename InclusionIndicator >
			double compute_sum_of_variances( Eigen::VectorXd const& levels, Probabilities const& probabilities, Eigen::VectorXd const& fallback, InclusionIndicator const& inclusion ) {
				assert( levels.size() == fallback.size() ) ;
				assert( levels.size() == probabilities.cols() ) ;
				assert( inclusion.size() == probabilities.rows() ) ;
				Eigen::VectorXd levels_squared = ( levels.array() * levels.array() ) ;

				double result = 0.0 ;
				for( int i = 0; i < probabilities.rows(); ++i ) {
					if( inclusion( i ) == 1 ) {
						double const c = probabilities.row( i ).sum() ;
						result += compute_variance( levels, levels_squared, probabilities.row( i ).transpose() + ( 1 - c ) * fallback ) ;
					}
				}
				return result ;
			}
	
		}
		
		InfoComputation::InfoComputation():
			m_diploid_fallback_distribution( Eigen::VectorXd::Zero( 3 )),
			m_haploid_fallback_distribution( Eigen::VectorXd::Zero( 2 )),
			m_diploid_levels( Eigen::VectorXd::LinSpaced( 3, 0, 2 )),
			m_haploid_levels( Eigen::VectorXd::LinSpaced( 2, 0, 1 ))
		{}
		
		void InfoComputation::compute(
			VariantIdentifyingData const& snp,
			Genotypes const& genotypes,
			Ploidy const& ploidy
		) {
			compute(
				snp, genotypes, ploidy,
				std::vector< metro::SampleRange >( 1, metro::SampleRange( 0, genotypes.rows()) )
			) ;
		}
			
		void InfoComputation::compute(
			VariantIdentifyingData const& snp,
			Genotypes const& genotypes,
			Ploidy const& ploidy,
			std::vector< metro::SampleRange > const& included_samples
		) {
			m_info = std::numeric_limits< double >::quiet_NaN() ;
			m_impute_info = std::numeric_limits< double >::quiet_NaN() ;

			// we don't compute for multiallelics currently
			if( snp.number_of_alleles() == 2 ) {
				compute_impl( snp, genotypes, ploidy, included_samples ) ;
			}
		}
	
		void InfoComputation::compute_impl(
			VariantIdentifyingData const& snp,
			Genotypes const& genotypes,
			Ploidy const& ploidy,
			std::vector< metro::SampleRange > const& included_samples
		) {
			typedef Eigen::Block< Genotypes const > GenotypeBlock ;
			typedef Eigen::VectorBlock< Eigen::VectorXd const > PloidyBlock ;

			// This function assumes every sample is haploid or diploid.
			Eigen::VectorXd const haploids = ( ploidy.array() == 1 ).cast< double >() ;
			Eigen::VectorXd const diploids = ( ploidy.array() == 2 ).cast< double >() ;
			
			double b_allele_count = 0 ;
			double total_allele_count = 0 ;
			double autosomal_b_allele_count = 0 ;		// for computations that pretend ploidy == 2 everywhere.
			double autosomal_total_allele_count = 0 ;	// for computations that pretend ploidy == 2 everywhere.

			double number_of_haploid_info_samples = 0 ;
			double number_of_diploid_info_samples = 0 ;
			double total_probability = 0 ;
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				GenotypeBlock const block = genotypes.block( included_samples[i].begin(), 0, included_samples[i].end() - included_samples[i].begin(), genotypes.cols() ) ;
				PloidyBlock const haploid_block = haploids.segment( included_samples[i].begin(), included_samples[i].size() ) ;
				PloidyBlock const diploid_block = diploids.segment( included_samples[i].begin(), included_samples[i].size() ) ;
		
				b_allele_count += (( block.col( 1 ) + 2.0 * block.col( 2 ) ).array() * diploid_block.array() ).sum() ;
				b_allele_count += (( block.col( 1 ) ).array() * haploid_block.array() ).sum() ;

				total_allele_count += 2 * ( block.rowwise().sum().array() * diploid_block.array() ).sum() ;
				total_allele_count += ( block.rowwise().sum().array() * haploid_block.array() ).sum() ;

				autosomal_b_allele_count += ( block.col( 1 ) + 2.0 * block.col( 2 ) ).array().sum() ;
				autosomal_total_allele_count += 2 * block.sum() ;
			
				number_of_haploid_info_samples += haploid_block.sum() ;
				number_of_diploid_info_samples += diploid_block.sum() ;
				total_probability += block.sum() ;
			}
			
			// MLE estimate of allele frequency
			double const theta_mle = b_allele_count / total_allele_count ;
			double const autosomal_theta_mle = autosomal_b_allele_count / autosomal_total_allele_count ;
			double const theta_est = theta_mle ;

			m_diploid_fallback_distribution( 0 ) = ( 1 - theta_mle ) * ( 1 - theta_mle ) ;
			m_diploid_fallback_distribution( 1 ) = 2.0 * theta_mle * ( 1 - theta_mle ) ;
			m_diploid_fallback_distribution( 2 ) = theta_mle * theta_mle ;
		
			m_haploid_fallback_distribution( 0 ) = 1 - theta_mle ;
			m_haploid_fallback_distribution( 1 ) = theta_mle ;
		
			//std::cerr << "theta = " << theta_mle << ", fallback_distribution = " << fallback_distribution.transpose() << ".\n" ;
		
			double diploid_info_term = 0.0 ;
			double haploid_info_term = 0.0 ;
			double impute_info_term = 0.0 ;
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				GenotypeBlock const diploid_genotype_block = genotypes.block( included_samples[i].begin(), 0, included_samples[i].size(), genotypes.cols() ) ;
				GenotypeBlock const haploid_genotype_block = genotypes.block( included_samples[i].begin(), 0, included_samples[i].size(), 2 ) ;
				PloidyBlock const haploid_block = haploids.segment( included_samples[i].begin(), included_samples[i].size() ) ;
				PloidyBlock const diploid_block = diploids.segment( included_samples[i].begin(), included_samples[i].size() ) ;

				haploid_info_term -= compute_sum_of_variances(
					m_haploid_levels, haploid_genotype_block,
					m_haploid_fallback_distribution, haploid_block
				)
				/ ( theta_mle * ( 1 - theta_mle ) ) ;

				diploid_info_term -= compute_sum_of_variances(
					m_diploid_levels, diploid_genotype_block,
					m_diploid_fallback_distribution, diploid_block
				)
				/ ( 2.0 * theta_mle * ( 1 - theta_mle ) ) ;
			
				impute_info_term -= compute_sum_of_variances(
					m_diploid_levels, diploid_genotype_block,
					Eigen::VectorXd::Zero( 3 ), Eigen::VectorXd::Constant( diploid_genotype_block.rows(), 1 )
				)
				/ ( 2.0 * autosomal_theta_mle * ( 1 - autosomal_theta_mle ) ) ;
			}
		
			double info = std::numeric_limits< double >::quiet_NaN() ;
			double impute_info = std::numeric_limits< double >::quiet_NaN() ;

			if(( number_of_haploid_info_samples + number_of_diploid_info_samples ) > 0 ) {
				info = 1.0 + (haploid_info_term + diploid_info_term) / ( number_of_haploid_info_samples + number_of_diploid_info_samples ) ;
			}
			
			if( total_probability > 0 ) {
				impute_info = 1.0 + impute_info_term / total_probability ;
			}

			m_info = info ;
			m_impute_info = impute_info ;
		}

	}
}
