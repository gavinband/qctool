#include <string>
#include <boost/signals2/signal.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "components/CallComparerComponent/PairwiseCallComparer.hpp"
#include "components/CallComparerComponent/PairwiseCallComparerManager.hpp"

PairwiseCallComparerManager::UniquePtr PairwiseCallComparerManager::create() {
	PairwiseCallComparerManager::UniquePtr result( new PairwiseCallComparerManager() ) ;
	return result ;
}

PairwiseCallComparerManager::PairwiseCallComparerManager()
{}

void PairwiseCallComparerManager::add_comparer( std::string const& name, PairwiseCallComparer::UniquePtr comparer ) {
	m_comparers.insert( name, comparer ) ;
}

void PairwiseCallComparerManager::send_comparisons_to( ComparisonClient::SharedPtr client ) {
	m_comparison_clients.push_back( client ) ;
}

void PairwiseCallComparerManager::send_merge_to( MergeClient::SharedPtr client ) {
	m_merge_clients.push_back( client ) ;
}

void PairwiseCallComparerManager::set_merger( Merger::UniquePtr merger ) {
	m_merger = merger ;
}

void PairwiseCallComparerManager::begin_processing_snp( genfile::SNPIdentifyingData const& snp ) {
	m_calls.clear() ;
	m_snp = snp ;
}

void PairwiseCallComparerManager::add_calls( std::string const& name, genfile::SingleSNPGenotypeProbabilities const& calls ) {
	m_calls[ name ] = calls ;
}

void PairwiseCallComparerManager::end_processing_snp() {
	for( std::size_t i = 0; i < m_comparison_clients.size(); ++i ) {
		m_comparison_clients[i]->begin_comparisons( m_snp ) ;
	}
	if( m_merger.get() ) {
		m_merger->begin_comparisons( m_snp ) ;
	}
	
	for( Calls::const_iterator i = m_calls.begin(); i != m_calls.end(); ++i ) {
		Calls::const_iterator j = i ;
		for( ++j; j != m_calls.end(); ++j ) {
			for( Comparers::const_iterator comparer_i = m_comparers.begin(); comparer_i != m_comparers.end(); ++comparer_i ) {
				std::map< std::string, genfile::VariantEntry > const& result = comparer_i->second->compare( i->second, j->second ) ;
				std::map< std::string, genfile::VariantEntry >::const_iterator
					result_i = result.begin(),
					result_end = result.end() ;
				for( ; result_i != result_end; ++result_i ) {
					send_comparisons_to_clients( i->first, j->first, comparer_i->first, result_i->first, result_i->second ) ;
					if( m_merger.get() ) {
						m_merger->set_result( i->first, j->first, comparer_i->first, result_i->first, result_i->second ) ;
					}
				}
			}
		}
	}
	
	for( std::size_t i = 0; i < m_comparison_clients.size(); ++i ) {
		m_comparison_clients[i]->end_comparisons() ;
	}

	if( m_merger.get() ) {
		m_merger->end_comparisons() ;
		for( std::size_t i = 0; i < m_merge_clients.size(); ++i ) {
			m_merge_clients[i]->begin_comparisons( m_snp ) ;
		}
		send_merge_to_clients( m_merger->get_spec(), "concordant_calls", m_merger->get_result_as_string() ) ;
		for( std::size_t i = 0; i < m_merge_clients.size(); ++i ) {
			m_merge_clients[i]->end_comparisons() ;
		}
	}
}

void PairwiseCallComparerManager::send_comparisons_to_clients(
	std::string const& first_callset,
	std::string const& second_callset,
	std::string const& comparison,
	std::string const& variable,
	genfile::VariantEntry const& value
) {
	for( std::size_t i = 0; i < m_comparison_clients.size(); ++i ) {
		m_comparison_clients[i]->set_result( first_callset, second_callset, comparison, variable, value ) ;
	}
}

void PairwiseCallComparerManager::send_merge_to_clients(
	std::string const& comparison,
	std::string const& variable,
	genfile::VariantEntry const& value
) {
	for( std::size_t i = 0; i < m_merge_clients.size(); ++i ) {
		m_merge_clients[i]->set_result( comparison, variable, value ) ;
	}
}
