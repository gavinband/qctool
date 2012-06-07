
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
#include "components/CallComparerComponent/LeastMissingConsensusCaller.hpp"

LeastMissingConsensusCaller::LeastMissingConsensusCaller():
	m_call_threshhold( 0.9 )
{}

void LeastMissingConsensusCaller::begin_processing_snps( std::size_t number_of_samples ) {
	m_number_of_samples = number_of_samples ;
}

void LeastMissingConsensusCaller::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	std::vector< std::string > const& call_names = get_consensus_call_names() ;
	m_genotypes.resize( call_names.size() ) ;
	std::map< std::string, std::vector< genfile::VariantEntry > > info ;
	std::size_t chosen_call = 0 ;
	double chosen_call_missingness = std::numeric_limits< double >::max() ;
	if( call_names.size() > 0 ) {
		for( std::size_t i = 0; i < call_names.size(); ++i ) {
			{
				genfile::vcf::GenotypeSetter< Eigen::MatrixBase< Eigen::MatrixXd > > setter( m_genotypes[i] ) ;
				data_reader.get( call_names[i], setter ) ;
			}
			double missing_calls = 0 ;
			for( int j = 0; j < m_genotypes[i].rows(); ++j ) {
				if( m_genotypes[i].row( j ).maxCoeff() < m_call_threshhold ) {
					++missing_calls ;
				}
			}
			info[ call_names[i] ].push_back( missing_calls / m_genotypes[i].rows() ) ;
			if( missing_calls < chosen_call_missingness ) {
				chosen_call = i ;
				chosen_call_missingness = missing_calls ;
			}
		}
	
		info[ "call" ].push_back( call_names[ chosen_call ] ) ;

		send_results(
			snp,
			m_genotypes[ chosen_call ],
			info
		) ;
	} else {
		m_genotypes.resize( 1 ) ;
		m_genotypes[0].setZero( m_number_of_samples, 3 ) ;
		info[ "call" ].push_back( genfile::MissingValue() ) ;
		send_results(
			snp,
			m_genotypes[ 0 ],
			info
		) ;
	}
}

