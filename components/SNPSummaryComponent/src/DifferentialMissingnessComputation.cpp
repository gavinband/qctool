#include <cassert>
#include <boost/math/distributions/chi_squared.hpp>
#include <Eigen/Core>
#include "genfile/Error.hpp"
#include "components/SNPSummaryComponent/DifferentialMissingnessComputation.hpp"
#include "metro/FishersExactTest.hpp"
#include "metro/likelihood/Multinomial.hpp"
#include "metro/likelihood/ProductOfMultinomials.hpp"

DifferentialMissingnessComputation::UniquePtr DifferentialMissingnessComputation::create( std::string const& stratification_name, StrataMembers const& strata_members ) {
	return DifferentialMissingnessComputation::UniquePtr(
		new DifferentialMissingnessComputation( stratification_name, strata_members )
	) ;
}

DifferentialMissingnessComputation::DifferentialMissingnessComputation( std::string const& stratification_name, StrataMembers const& strata_members, double threshhold ):
	m_stratification_name( stratification_name ),
	m_strata_members( strata_members ),
	m_strata_levels( compute_strata_levels( m_strata_members ) ),
	m_threshhold( threshhold )
{
	if( !strata_members.size() == 2 ) {
		throw genfile::BadArgumentError( "DifferentialMissingnessComputation::DifferentialMissingnessComputation()", "strata_members( size " + genfile::string_utils::to_string( strata_members.size() ) + ")" ) ;
	}
}

std::vector< int > DifferentialMissingnessComputation::compute_strata_levels( StrataMembers const& strata_members ) const {
	std::vector< int > result ;
	int level = 0 ;
	for( StrataMembers::const_iterator i = strata_members.begin(); i != strata_members.end(); ++i, ++level ) {
		for( std::size_t j = 0; j < i->second.size(); ++j ) {
			std::size_t const index = i->second[j] ;
			result.resize( std::max( result.size(), index+1 ), -1 ) ;
			result[ index ] = level ;
		}
	}
	return result ;
}

void DifferentialMissingnessComputation::operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, SampleSexes const&, genfile::VariantDataReader&, ResultCallback callback ) {
	// construct a table
	// 
	//                 missing     not missing
	//   stratum 1       a              b
	//   stratum 2       c              d
	// 
	// on which we can do a Fisher's exact test or a product of binomial-type test.
	//
	int const number_of_strata = 1 + *std::max_element( m_strata_levels.begin(), m_strata_levels.end() ) ;
	Eigen::MatrixXd table( number_of_strata, 2 ) ;
	table.setZero() ;
	for( int i = 0; i < genotypes.rows(); ++i ) {
		if( m_strata_levels[i] >= 0 ) {
			if( genotypes.row( i ).sum() < m_threshhold ) {
				table( m_strata_levels[i], 0 )++ ;
			}
			else {
				table( m_strata_levels[i], 1 )++ ;
			}		
		}
	}
	
	{
		StrataMembers::const_iterator i = m_strata_members.begin() ;
		StrataMembers::const_iterator const end_i = m_strata_members.end() ;
		for( int level = 0; i != end_i; ++level, ++i ) {
			callback( "missing_call_proportion [" + m_stratification_name + "=" + genfile::string_utils::to_string( i->first ) + "]", table( level, 0 ) / table.row( level ).sum() )  ;
		}

	}
	
	if( table.row(0).sum() > 0 && table.row(1).sum() > 0 ) {
		std::string const stub = "missingness_by_" + m_stratification_name ;
		try {
			metro::FishersExactTest test( table ) ;
			callback( stub + "_exact_pvalue", test.get_pvalue() )  ;
			callback( stub + "_sample_odds_ratio", test.get_OR() )  ;
		}
		catch( std::exception const& e ) {
			//
		}
		
		{
			metro::likelihood::Multinomial< double, Eigen::VectorXd, Eigen::MatrixXd > null_model( table.colwise().sum() ) ;
			null_model.evaluate_at( null_model.get_MLE() ) ;
			metro::likelihood::ProductOfMultinomials< double, Eigen::VectorXd, Eigen::MatrixXd > alternative_model( table ) ;
			alternative_model.evaluate_at( alternative_model.get_MLE() ) ;
			double likelihood_ratio_statistic = -2.0 * ( null_model.get_value_of_function() - alternative_model.get_value_of_function() ) ;

			if( likelihood_ratio_statistic != likelihood_ratio_statistic || likelihood_ratio_statistic < 0.0 ) {
				likelihood_ratio_statistic = std::numeric_limits< double >::quiet_NaN() ;
			}
			else {
				boost::math::chi_squared_distribution< double > chi_squared( table.rows() - 1 ) ;

				double p_value = boost::math::cdf(
					boost::math::complement(
						chi_squared,
						likelihood_ratio_statistic
					)
				) ;
				
				
				callback( stub + "_lr_pvalue", p_value ) ;
				callback( stub + "_lr_statistic", likelihood_ratio_statistic ) ;
			}
		}
	}
}

std::string DifferentialMissingnessComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
	return prefix + "DifferentialMissingnessComputation" ;
}
