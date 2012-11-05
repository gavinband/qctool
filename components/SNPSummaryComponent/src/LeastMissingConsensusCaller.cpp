
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>
#include <map>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "components/SNPSummaryComponent/LeastMissingConsensusCaller.hpp"

LeastMissingConsensusCaller::LeastMissingConsensusCaller():
	m_call_threshhold( 0.9 )
{}

void LeastMissingConsensusCaller::set_result(
	std::string const& comparison,
	std::string const& accepted_calls,
	PairwiseCallComparerManager::Calls const& calls
) {
	std::vector< std::string > const& call_names = genfile::string_utils::split( accepted_calls, "," ) ;
	m_genotypes.resize( call_names.size() ) ;
	std::map< std::string, std::vector< genfile::VariantEntry > > info ;
	std::string chosen_call = "" ;
	double chosen_call_missingness = std::numeric_limits< double >::max() ;
	if( call_names.size() > 0 ) {
		for( std::size_t i = 0; i < call_names.size(); ++i ) {
			Eigen::MatrixXd const& genotype_calls = calls.find( call_names[i] )->second ;
			double missing_calls = 0 ;
			for( int j = 0; j < genotype_calls.rows(); ++j ) {
				if( genotype_calls.row( j ).maxCoeff() < m_call_threshhold ) {
					++missing_calls ;
				}
			}
			info[ call_names[i] ].push_back( missing_calls / genotype_calls.rows() ) ;
			if( missing_calls < chosen_call_missingness ) {
				chosen_call = call_names[i] ;
				chosen_call_missingness = missing_calls ;
			}
		}
	
		info[ "call" ].push_back( chosen_call ) ;

		send_results(
			get_snp(),
			calls.find( chosen_call )->second,
			info
		) ;
	} else {
		info[ "call" ].push_back( genfile::MissingValue() ) ;
		send_results(
			get_snp(),
			Eigen::MatrixXd::Zero( get_number_of_samples(), 3 ),
			info
		) ;
	}
}

