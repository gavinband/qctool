
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include "genfile/Error.hpp"
#include "components/SNPSummaryComponent/SNPHWE.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"

namespace {
	struct AlleleProportionComputation: public SNPSummaryComputation
	{
		void operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, genfile::VariantDataReader&, ResultCallback callback ) {
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
		
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const { return prefix + "AlleleProportionComputation" ; }
	} ;
	
	struct HWEComputation: public SNPSummaryComputation
	{
		void operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, genfile::VariantDataReader&, ResultCallback callback ) {
			double const
				AA = std::floor( genotypes.col(0).sum() + 0.5 ),
				AB = std::floor( genotypes.col(1).sum() + 0.5 ),
				BB = std::floor( genotypes.col(2).sum() + 0.5 ) ;

			if( AA > 0.5 || AB > 0.5 || BB > 0.5 ) {
				double HWE_pvalue = SNPHWE( AB, AA, BB ) ;
				callback( "minus_log10_exact_HW_p_value", -log10( HWE_pvalue ) ) ;
			}
			else {
				callback( "minus_log10_exact_HW_p_value", genfile::MissingValue() ) ;
			}
		}
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const { return prefix + "HWEComputation" ; }
	} ;
	
	struct MissingnessComputation: public SNPSummaryComputation {
		MissingnessComputation( double call_threshhold = 0.9 ): m_call_threshhold( call_threshhold ) {}
		void operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, genfile::VariantDataReader&, ResultCallback callback ) {
			double missingness = double( genotypes.rows() ) - genotypes.array().sum() ;
			callback( "missing proportion", missingness / double( genotypes.rows() ) ) ;
			
			double missing_calls = 0.0 ;
			for( int i = 0; i < genotypes.rows(); ++i ) {
				if( genotypes.row( i ).maxCoeff() < m_call_threshhold ) {
					++missing_calls ;
				}
			}
			callback( "missing call proportion", missing_calls / double( genotypes.rows() )) ;

			callback( "AA", genotypes.col(0).sum() ) ;
			callback( "AB", genotypes.col(1).sum() ) ;
			callback( "BB", genotypes.col(2).sum() ) ;
		}
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const { return prefix + "MissingnessComputation" ; }
	private:
		double const m_call_threshhold ;
	} ;

	struct InformationComputation: public SNPSummaryComputation {
		void operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, genfile::VariantDataReader&, ResultCallback callback ) {
			Eigen::VectorXd const e = genotypes.col(1) + (2.0 * genotypes.col(2) ) ;
			Eigen::VectorXd const f = genotypes.col(1) + (4.0 * genotypes.col(2) ) ;
			double const non_missingness = genotypes.sum() ;
			Eigen::VectorXd adjustment1 = Eigen::VectorXd::Ones( genotypes.rows() ) - genotypes.rowwise().sum() ;
			Eigen::VectorXd const adjustment2 = adjustment1.array().square() ;
			adjustment1.array() *= e.array() ;

			if( non_missingness == 0.0 ) {
				callback( "info", genfile::MissingValue() ) ;
				callback( "impute_info", genfile::MissingValue() ) ;
			}
			else {
				double const theta_mle = e.sum() / ( 2.0 * non_missingness ) ;

				if( theta_mle == 0.0 || theta_mle == 1.0 ) {
					callback( "info", 1.0 ) ;
					callback( "impute_info", 1.0 ) ;
				}
				else {
					double const variance = ( f.array() - ( e.array().square() ) ).sum() ;
					double adjustment = 4.0 * theta_mle * adjustment1.sum() + 4.0 * theta_mle * theta_mle * adjustment2.sum() ;

					adjustment += ( genotypes.rows() - non_missingness ) * 2.0 * theta_mle * ( 1.0 + theta_mle ) ;

					double denominator = 2.0 * genotypes.rows() * theta_mle * ( 1.0 - theta_mle ) ;

					callback( "info", 1.0 - ( ( variance + adjustment ) / denominator ) ) ;
					callback( "impute_info", 1.0 - ( variance / denominator ) ) ;
				}
			}
		}
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const { return prefix + "InformationComputation" ; }
	} ;

}

SNPSummaryComputation::UniquePtr SNPSummaryComputation::create(
	std::string const& name
) {
	UniquePtr result ;
	if( name == "alleles" ) { result.reset( new AlleleProportionComputation()) ; }
	else if( name == "HWE" ) { result.reset( new HWEComputation()) ; }
	else if( name == "missingness" ) { result.reset( new MissingnessComputation()) ; }
	else if( name == "information" ) { result.reset( new InformationComputation()) ; }
	else {
		throw genfile::BadArgumentError( "SNPSummaryComputation::create()", "name=\"" + name + "\"" ) ;
	}
	return result ;
}
