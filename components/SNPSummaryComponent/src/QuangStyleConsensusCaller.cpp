
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
#include "components/SNPSummaryComponent/QuangStyleConsensusCaller.hpp"

// #define DEBUG_QUANGSTYLECONSENSUSCALLER 1
QuangStyleConsensusCaller::QuangStyleConsensusCaller( double threshhold, unsigned int minimum_consensus ):
	m_call_threshhold( threshhold ),
	m_minimum_consensus( minimum_consensus )
{}

void QuangStyleConsensusCaller::set_result(
	std::string const& comparison,
	std::string const& accepted_calls,
	PairwiseCallComparerManager::Calls const& calls
) {
	std::vector< std::string > call_names = genfile::string_utils::split( accepted_calls, "," ) ;
	m_genotypes.resize( call_names.size() ) ;
	std::map< std::string, std::vector< genfile::VariantEntry > > info ;

	if( call_names.size() > 0 ) {
		m_consensus_counts.setConstant( get_number_of_samples(), call_names.size() ) ;
		m_result_calls.setZero( get_number_of_samples() ) ;

		for( std::size_t i = 0; i < call_names.size(); ++i ) {
#if DEBUG_QUANGSTYLECONSENSUSCALLER
			std::cerr << "QuangStyleConsensusCaller::processed_snp(): looking at call " << call_names[i] << "...\n" ;
#endif
			// turn genotype probs into hard-called genotypes, 0 = missing, 1, 2, 3 as levels.
			{
				Eigen::MatrixXd const& genotype_probs = calls.find( call_names[i] )->second ;

				assert( genotype_probs.rows() == get_number_of_samples() ) ;
				m_genotypes.resize( genotype_probs.rows() ) ;
				{
					genfile::vcf::impl::ThreshholdedCallGetter1< Eigen::VectorXd > setter = genfile::vcf::get_threshholded_calls( m_genotypes, m_call_threshhold, 0, 1, 2, 3 ) ;
					setter.set_number_of_samples( get_number_of_samples() ) ;
					for( std::size_t sample = 0; sample < get_number_of_samples(); ++sample ) {
						setter.set( sample, genotype_probs( sample, 0 ), genotype_probs( sample, 1 ), genotype_probs( sample, 2 )) ;
					}
				}
				assert( std::size_t( m_genotypes.size() ) == get_number_of_samples() ) ;
			}
			
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
		m_result_calls.array() *= ( m_consensus_counts.array() > m_minimum_consensus ).cast< double >() ;

		info[ "conflicting_calls" ].push_back( ( m_result_calls.array() < 0 ).cast< int >().sum() ) ;
		info[ "missing_calls" ].push_back( ( m_result_calls.array() == 0.0 ).cast< int >().sum() ) ;

		// Convert to probabilities.
		m_result_probs.setZero( get_number_of_samples(), 3 ) ;
		for( std::size_t i = 0; i < get_number_of_samples(); ++i ) {
			if( m_result_calls( i ) > 0.0 ) {
				m_result_probs( i, m_result_calls( i ) - 1 ) = 1 ;
			}
		}

		send_results( get_snp(), m_result_probs, info ) ;
	} else {
		m_result_probs.setZero( get_number_of_samples(), 3 ) ;
		send_results( get_snp(), m_result_probs, info ) ;
	}
}

