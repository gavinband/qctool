#include <string>
#include <vector>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include "genfile/Error.hpp"
#include "SNPHWE.hpp"
#include "SNPSummaryComputation.hpp"

namespace {
	struct AlleleProportionComputation: public SNPSummaryComputation
	{
		void operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, ResultCallback callback ) const {
			double result = std::numeric_limits< double >::quiet_NaN() ;
			double const a_allele_freq = ( ( 2.0 * genotypes.col(0).sum() ) + genotypes.col(1).sum() ) / ( 2.0 * genotypes.sum() ) ;
			double const b_allele_freq = ( ( 2.0 * genotypes.col(2).sum() ) + genotypes.col(1).sum() ) / ( 2.0 * genotypes.sum() ) ;

			bool minor_is_first = a_allele_freq <= b_allele_freq ;

			callback( "alleleA_frequency", a_allele_freq ) ;
			callback( "alleleB_frequency", b_allele_freq ) ;

			if( a_allele_freq < b_allele_freq ) {
				callback( "minor_allele_frequency", a_allele_freq ) ;
				callback( "minor_allele", snp.get_first_allele() ) ;
			}
			else if( a_allele_freq > b_allele_freq ) {
				callback( "minor_allele_frequency", b_allele_freq ) ;
				callback( "minor_allele", snp.get_second_allele() ) ;
			}
		}
	} ;
	
	struct HWEComputation: public SNPSummaryComputation
	{
		void operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, ResultCallback callback ) const {
			double const
				AA = std::floor( genotypes.col(0).sum() + 0.5 ),
				AB = std::floor( genotypes.col(1).sum() + 0.5 ),
				BB = std::floor( genotypes.col(2).sum() + 0.5 ) ;

			if( AA > 0.5 || AB > 0.5 || BB > 0.5 ) {
				double HWE_pvalue = SNPHWE( AB, AA, BB ) ;
				callback( "-log_10_exact_HW_p-value", HWE_pvalue ) ;
			}
			else {
				callback( "-log_10_exact_HW_p-value", genfile::MissingValue() ) ;
			}
		}
	} ;
	
	struct MissingnessComputation: public SNPSummaryComputation {
		MissingnessComputation( double call_threshhold = 0.9 ): m_call_threshhold( call_threshhold ) {}
		void operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, ResultCallback callback ) const {
			double missingness = genotypes.rows() - genotypes.sum() ;
			callback( "missing proportion", missingness / double( genotypes.rows() ) ) ;
			callback( "missingness", missingness ) ;
			callback( "non-missingness", genotypes.rows() - missingness ) ;
			
			double missing_calls = 0.0 ;
			for( int i = 0; i < genotypes.rows(); ++i ) {
				if( genotypes.row( i ).maxCoeff() < m_call_threshhold ) {
					++missing_calls ;
				}
			}
			callback( "missing calls", missing_calls ) ;
			
		}
	private:
		double const m_call_threshhold ;
	} ;

	struct InformationComputation: public SNPSummaryComputation {
		void operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, ResultCallback callback ) const {
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

				double const variance = ( f.array() - ( e.array().square() ) ).sum() ;
				double adjustment = 4.0 * theta_mle * adjustment1.sum() + 4.0 * theta_mle * theta_mle * adjustment2.sum() ;

				double const missingness = genotypes.rows() - non_missingness ;
				adjustment += ( genotypes.rows() - non_missingness ) * 2.0 * theta_mle * ( 1.0 + theta_mle ) ;

				std::cerr << snp.get_SNPID() << ": "
					<< "genotypes: " << genotypes.rows() << "x" << genotypes.cols() << ", "
					<< "sum: " << genotypes.sum() << ", "
					<< "e.sum() == " << e.sum()  << ", "
					<< "theta_mle = " << theta_mle << ", "
					<< "variance = " << variance << ", "
					<< "adjustment = " << adjustment << ", "
					<< ".\n" ;

				double denominator = 2.0 * genotypes.rows() * theta_mle * ( 1.0 - theta_mle ) ;

				callback( "info", 1.0 - ( ( variance + adjustment ) / denominator ) ) ;
				callback( "impute_info", 1.0 - ( variance / denominator ) ) ;
			}
		}
	} ;

}

void SNPSummaryComputation::list_computations( boost::function< void ( std::string ) > callback ) {
	callback( "SNPID" ) ;
	callback( "rsid" ) ;
	callback( "alleles" ) ;
	callback( "HWE" ) ;
	callback( "missingness" ) ;
	callback( "information" ) ;
}


SNPSummaryComputation::UniquePtr SNPSummaryComputation::create( std::string const& name ) {
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
