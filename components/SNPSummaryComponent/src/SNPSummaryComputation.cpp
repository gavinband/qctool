
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
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "components/SNPSummaryComponent/HWEComputation.hpp"

namespace snp_summary_component {
	struct AlleleFrequencyComputation: public SNPSummaryComputation
	{
		void operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, SampleSexes const& sexes, genfile::VariantDataReader&, ResultCallback callback ) {
			genfile::Chromosome const& chromosome = snp.get_position().chromosome() ;
			if( chromosome.is_sex_determining() ) {
				compute_sex_chromosome_frequency( snp, genotypes, sexes, callback ) ;
			}
			else {
				compute_autosomal_frequency( snp, genotypes, callback ) ;
			}
		}
		
		void compute_sex_chromosome_frequency( SNPIdentifyingData const& snp, Genotypes const& genotypes, SampleSexes const& sexes, ResultCallback callback ) {
			Genotypes male_genotypes = genotypes ;
			Genotypes female_genotypes = genotypes ;
			assert( std::size_t( genotypes.rows() ) == sexes.size() ) ;

			for( int i = 0; i < genotypes.rows(); ++i ) {
				if( sexes[i] != 'm' ) {
					male_genotypes.row(i).setZero() ;
				} else if( sexes[i] != 'f' ) {
					female_genotypes.row(i).setZero() ;
				}
			}
			
			double const b_allele_count = male_genotypes.col(1).sum()
				+ ( ( 2.0 * female_genotypes.col(2).sum() ) + female_genotypes.col(1).sum() ) ;

			double const b_allele_freq = b_allele_count / ( male_genotypes.sum() + 2.0 * female_genotypes.sum() ) ;

			double const a_allele_freq = 1 - b_allele_freq ;

			callback( "alleleA_frequency", a_allele_freq ) ;
			callback( "alleleB_frequency", b_allele_freq ) ;

			if( a_allele_freq < b_allele_freq ) {
				callback( "minor_allele_frequency", a_allele_freq ) ;
				callback( "minor_allele", snp.get_first_allele() ) ;
				callback( "major_allele", snp.get_second_allele() ) ;
			}
			else if( a_allele_freq > b_allele_freq ) {
				callback( "minor_allele_frequency", b_allele_freq ) ;
				callback( "minor_allele", snp.get_second_allele() ) ;
				callback( "major_allele", snp.get_first_allele() ) ;
			} else {
				callback( "minor_allele_frequency", a_allele_freq ) ;
			}
		}

		void compute_autosomal_frequency( SNPIdentifyingData const& snp, Genotypes const& genotypes, ResultCallback callback ) {
			double const a_allele_freq = ( ( 2.0 * genotypes.col(0).sum() ) + genotypes.col(1).sum() ) / ( 2.0 * genotypes.sum() ) ;
			double const b_allele_freq = ( ( 2.0 * genotypes.col(2).sum() ) + genotypes.col(1).sum() ) / ( 2.0 * genotypes.sum() ) ;

			callback( "alleleA_frequency", a_allele_freq ) ;
			callback( "alleleB_frequency", b_allele_freq ) ;

			if( a_allele_freq < b_allele_freq ) {
				callback( "minor_allele_frequency", a_allele_freq ) ;
				callback( "minor_allele", snp.get_first_allele() ) ;
				callback( "major_allele", snp.get_second_allele() ) ;
			}
			else if( a_allele_freq > b_allele_freq ) {
				callback( "minor_allele_frequency", b_allele_freq ) ;
				callback( "minor_allele", snp.get_second_allele() ) ;
				callback( "major_allele", snp.get_first_allele() ) ;
			} else {
				callback( "minor_allele_frequency", a_allele_freq ) ;
			}
		}
		
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const { return prefix + "AlleleFrequencyComputation" ; }
	} ;
	
	struct MissingnessComputation: public SNPSummaryComputation {
		MissingnessComputation( double call_threshhold = 0.9 ): m_call_threshhold( call_threshhold ) {}
		void operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, SampleSexes const& sexes, genfile::VariantDataReader&, ResultCallback callback ) {
			double missingness = double( genotypes.rows() ) - genotypes.array().sum() ;
			callback( "missing proportion", missingness / double( genotypes.rows() ) ) ;
			
			double missing_calls = 0.0 ;
			for( int i = 0; i < genotypes.rows(); ++i ) {
				if( genotypes.row( i ).maxCoeff() < m_call_threshhold ) {
					++missing_calls ;
				}
			}
			callback( "missing call proportion", missing_calls / double( genotypes.rows() )) ;

			genfile::Chromosome const& chromosome = snp.get_position().chromosome() ;

			if( !chromosome.is_sex_determining() ) {
				callback( "AA", genotypes.col(0).sum() ) ;
				callback( "AB", genotypes.col(1).sum() ) ;
				callback( "BB", genotypes.col(2).sum() ) ;
			} else {
				compute_sex_chromosome_counts( snp, genotypes, sexes, callback ) ;
			}
		}
		
		void compute_sex_chromosome_counts( SNPIdentifyingData const& snp, Genotypes const& genotypes, SampleSexes const& sexes, ResultCallback callback ) {
			genfile::Chromosome const& chromosome = snp.get_position().chromosome() ;
			assert( chromosome == genfile::Chromosome( "0X" ) || chromosome == genfile::Chromosome( "0Y" )) ;
			assert( std::size_t( genotypes.rows() ) == sexes.size() ) ;

			std::map< char, Eigen::VectorXd > counts ;
			counts[ 'm' ] = Eigen::VectorXd::Zero( 3 ) ;
			counts[ 'f' ] = Eigen::VectorXd::Zero( 3 ) ;
			counts[ '.' ] = Eigen::VectorXd::Zero( 3 ) ;

			std::map< char, std::size_t > sample_counts ;

			for( std::size_t i = 0; i < sexes.size(); ++i ) {
				counts[ sexes[i] ] += genotypes.row( i ) ;
				++sample_counts[ sexes[i] ] ;
			}
			
			callback( "males_A", counts[ 'm' ]( 0 ) ) ;
			callback( "males_B", counts[ 'm' ]( 1 ) ) ;

			if( counts[ 'm' ]( 2 ) != 0 ) {
				callback( "males_BB", counts[ 'm' ]( 2 ) ) ;
				std::cerr << "!! ( MissingnessComputation::compute_sex_chromosome_counts() ): some males have BB probability!\n" ;
				throw genfile::BadArgumentError( " MissingnessComputation::compute_sex_chromosome_counts()", "genotypes" ) ;
			}

			if( chromosome == genfile::Chromosome( "0X" ) ) {
				callback( "females_AA", counts[ 'f' ]( 0 ) ) ;
				callback( "females_AB", counts[ 'f' ]( 1 ) ) ;
				callback( "females_BB", counts[ 'f' ]( 2 ) ) ;
			}
			
			if( sample_counts[ '.' ] > 0 ) {
				callback( "unknown_sex_AA", counts[ '.' ]( 0 ) ) ;
				callback( "unknown_sex_AB", counts[ '.' ]( 1 ) ) ;
				callback( "unknown_sex_BB", counts[ '.' ]( 2 ) ) ;
			}
		}
		
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const { return prefix + "MissingnessComputation" ; }
	private:
		double const m_call_threshhold ;
	} ;

	struct InformationComputation: public SNPSummaryComputation {
		void operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, SampleSexes const& sample_sexes, genfile::VariantDataReader&, ResultCallback callback ) {
			genfile::Chromosome const& chromosome = snp.get_position().chromosome() ;
			if( chromosome.is_sex_determining() ) {
				return ;
			}

			compute_info( snp, genotypes, sample_sexes, callback ) ;
		}
		
		void compute_info( SNPIdentifyingData const& snp, Genotypes const& genotypes, SampleSexes const&, ResultCallback callback ) {
			double theta_mle = ( genotypes.col( 1 ).sum() + 2.0 * genotypes.col( 2 ).sum() ) / ( 2.0 * genotypes.sum() ) ;
			
			Eigen::VectorXd const impute_fallback_distribution = Eigen::VectorXd::Zero( 3 ) ;
			Eigen::VectorXd fallback_distribution = Eigen::VectorXd::Zero( 3 ) ;
			fallback_distribution( 0 ) = ( 1 - theta_mle ) * ( 1 - theta_mle ) ;
			fallback_distribution( 1 ) = 2.0 * theta_mle * ( 1 - theta_mle ) ;
			fallback_distribution( 2 ) = theta_mle * theta_mle ;
			
			//std::cerr << "theta = " << theta_mle << ", fallback_distribution = " << fallback_distribution.transpose() << ".\n" ;
			
			Eigen::VectorXd const levels = Eigen::VectorXd::LinSpaced( 3, 0, 2 ) ;

			callback(
				"info",
				1.0 - (
					compute_expected_variance( levels, genotypes, fallback_distribution )
					/ ( 2.0 * theta_mle * ( 1 - theta_mle ) )
				)
			) ;
			callback(
				"impute_info",
				1.0 - (
					compute_expected_variance( levels, genotypes, impute_fallback_distribution )
					/ ( 2.0 * theta_mle * ( 1 - theta_mle ) )
				)
			) ;
		}

		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const { return prefix + "InformationComputation" ; }
		
	private:
		// treat the rows of the probabilities matrix as probabilities.
		// distribution for individual i is taken as a mixture of the distribution given by row i of probabilities,
		// and the fallback distribution if the sum of row i < 1, appropriately weighted.
		double compute_expected_variance( Eigen::VectorXd const& levels, Eigen::MatrixXd const& probabilities, Eigen::VectorXd const& fallback ) const {
			assert( levels.size() == fallback.size() ) ;
			assert( levels.size() == probabilities.cols() ) ;
			Eigen::VectorXd levels_squared = ( levels.array() * levels.array() ) ;

			double result = 0.0 ;
			for( int i = 0; i < probabilities.rows(); ++i ) {
				double const c = probabilities.row( i ).sum() ;
				result += compute_variance( levels, levels_squared, probabilities.row( i ).transpose() + ( 1 - c ) * fallback ) ;
			}
			return result / probabilities.rows() ;
		}
		
		double compute_variance( Eigen::VectorXd const& levels, Eigen::VectorXd const& levels_squared, Eigen::VectorXd const& probs ) const {
			double const mean = ( probs.transpose() * levels )(0) ;
			double const variance = ( probs.transpose() * levels_squared )(0) - ( mean * mean ) ;
			// std::cerr << "compute_variance: levels = " << levels.transpose() << ", levels_squared = " << levels_squared.transpose() << ", probs = " << probs.transpose() << ", variance = " << variance << ".\n" ;
			return variance ;
		}
		
	} ;

}

SNPSummaryComputation::UniquePtr SNPSummaryComputation::create(
	std::string const& name
) {
	UniquePtr result ;
	if( name == "alleles" ) { result.reset( new snp_summary_component::AlleleFrequencyComputation()) ; }
	else if( name == "HWE" ) { result.reset( new snp_summary_component::HWEComputation()) ; }
	else if( name == "missingness" ) { result.reset( new snp_summary_component::MissingnessComputation()) ; }
	else if( name == "information" ) { result.reset( new snp_summary_component::InformationComputation()) ; }
	else {
		throw genfile::BadArgumentError( "SNPSummaryComputation::create()", "name=\"" + name + "\"" ) ;
	}
	return result ;
}
