
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
#include "components/CallComparerComponent/QuangStyleConsensusCaller.hpp"

// #define DEBUG_QUANGSTYLECONSENSUSCALLER 1
QuangStyleConsensusCaller::QuangStyleConsensusCaller():
	m_call_threshhold( 0.9 )
{}

void QuangStyleConsensusCaller::begin_processing_snps( std::size_t number_of_samples ) {
	m_number_of_samples = number_of_samples ;
}

void QuangStyleConsensusCaller::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	std::vector< std::string > const& call_names = get_consensus_call_names() ;
	m_genotypes.resize( call_names.size() ) ;
	std::map< std::string, std::vector< genfile::VariantEntry > > info ;
	std::size_t chosen_call = 0 ;

	if( call_names.size() > 0 ) {
		m_consensus_counts.setConstant( m_number_of_samples, call_names.size() ) ;
		m_result_calls.setZero( m_number_of_samples ) ;

		for( std::size_t i = 0; i < call_names.size(); ++i ) {
#if DEBUG_QUANGSTYLECONSENSUSCALLER
			std::cerr << "QuangStyleConsensusCaller::processed_snp(): looking at call " << call_names[i] << "...\n" ;
#endif
			// get hard-called genotypes, 0 = missing, 1, 2, 3 as levels.
			{
				genfile::vcf::ThreshholdingGenotypeSetter< Eigen::VectorXd > setter( m_genotypes, m_call_threshhold, 0, 1, 2, 3 ) ;
				data_reader.get( call_names[i], setter ) ;
			}
			assert( std::size_t( m_genotypes.size() ) == m_number_of_samples ) ;

			// record the samples that are missing in this callset.
			m_consensus_counts.array() -= ( m_genotypes.array() == 0.0 ).cast< double >() ;

			if( i == 0 ) {
				m_result_calls = m_genotypes ;
			}
			else {
				// missing elts are encoded as 0, genotypes as 1, 2, 3.
				// we will encode conflict as large negative entries.
				// Apply the consensus call, i.e. set to 0 any that do not agree.
				// note that the term on the rhs might have 2s in it, but they only occur
				// if both the 0th and ith callsets are missing, when the result is 0 anyway.

				// replace any conflicting elements with something big and negative
				Eigen::VectorXd const conflicts =
					( m_result_calls.array() != 0.0 ).cast< double >()
					* ( m_genotypes.array() != 0.0 ).cast< double >()
					* ( m_genotypes.array() != m_result_calls.array() ).cast< double >() ;
				
				Eigen::VectorXd const good_calls = (
						( m_result_calls.array() == 0.0 ).cast< double >()
						* ( m_genotypes.array() != 0.0 ).cast< double >()
					) ;
				
				m_result_calls.array() += -10 * conflicts.array() ;
				m_result_calls.array() += good_calls.array() * m_genotypes.array() ;
			}
#if DEBUG_QUANGSTYLECONSENSUSCALLER
			std::cerr << "QuangStyleConsensusCaller::processed_snp(): result calls: " << m_result_calls.head( 10 ).transpose() << ".\n" ;
			std::cerr << "QuangStyleConsensusCaller::processed_snp(): genotypes: " << m_genotypes.head( 10 ).transpose() << ".\n" ;
			std::cerr << "QuangStyleConsensusCaller::processed_snp(): consensus_counts: " << m_consensus_counts.head( 10 ).transpose() << ".\n" ;
#endif
		}

		// set to missing anything with no consensus among algorithms.
		m_result_calls.array() *= ( m_consensus_counts.array() > 1.0 ).cast< double >() ;

		info[ "conflicting_calls" ].push_back( ( m_result_calls.array() < 0 ).cast< int >().sum() ) ;
		info[ "missing_calls" ].push_back( ( m_result_calls.array() == 0.0 ).cast< int >().sum() ) ;

		// Convert to probabilities.
		m_result_probs.setZero( m_number_of_samples, 3 ) ;
		for( std::size_t i = 0; i < m_number_of_samples; ++i ) {
			if( m_result_calls( i ) > 0.0 ) {
				m_result_probs( i, m_result_calls( i ) - 1 ) = 1 ;
			}
		}

		send_results( snp, m_result_probs, info ) ;
	} else {
		m_result_probs.setZero( m_number_of_samples, 3 ) ;
		send_results( snp, m_result_probs, info ) ;
	}
}

