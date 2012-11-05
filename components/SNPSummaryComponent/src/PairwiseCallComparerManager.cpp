
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <boost/signals2/signal.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "components/SNPSummaryComponent/PairwiseCallComparer.hpp"
#include "components/SNPSummaryComponent/PairwiseCallComparerManager.hpp"
#include "components/SNPSummaryComponent/FrequentistTestCallMerger.hpp"
#include "components/SNPSummaryComponent/AcceptAllCallMerger.hpp"

PairwiseCallComparerManager::UniquePtr PairwiseCallComparerManager::create() {
	PairwiseCallComparerManager::UniquePtr result( new PairwiseCallComparerManager() ) ;
	return result ;
}

PairwiseCallComparerManager::PairwiseCallComparerManager()
{}

void PairwiseCallComparerManager::add_comparer( std::string const& name, PairwiseCallComparer::UniquePtr comparer ) {
	if( comparer.get() ) {
		m_comparers.insert( name, comparer ) ;
	}
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

void PairwiseCallComparerManager::begin_processing_snps( std::size_t number_of_samples ) {
	for( std::size_t i = 0; i < m_merge_clients.size(); ++i ) {
		m_merge_clients[i]->begin_processing_snps( number_of_samples ) ;
	}
}

void PairwiseCallComparerManager::begin_processing_snp( genfile::SNPIdentifyingData const& snp ) {
	m_calls.clear() ;
	m_snp = snp ;
}

void PairwiseCallComparerManager::add_calls( std::string const& name, Eigen::MatrixXd const& calls ) {
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
		m_merger->add_callset( i->first ) ;
		Calls::const_iterator j = i ;
		for( ++j; j != m_calls.end(); ++j ) {
			for( Comparers::const_iterator comparer_i = m_comparers.begin(); comparer_i != m_comparers.end(); ++comparer_i ) {
				comparer_i->second->compare(
					i->second,
					j->second,
					boost::bind(
						&PairwiseCallComparerManager::send_comparisons_to_clients,
						this,
						i->first,
						j->first,
						comparer_i->first,
						_1,
						_2
					)
				) ;
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
		send_merge_to_clients( m_merger->get_spec(), m_merger->get_result_as_string(), m_calls ) ;
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
	if( m_merger.get() ) {
		m_merger->set_result( first_callset, second_callset, comparison, variable, value ) ;
	}
}

void PairwiseCallComparerManager::send_merge_to_clients(
	std::string const& comparison,
	std::string const& accepted_calls,
	Calls const& m_calls
) {
	for( std::size_t i = 0; i < m_merge_clients.size(); ++i ) {
		m_merge_clients[i]->set_result( comparison, accepted_calls, m_calls ) ;
	}
}

PairwiseCallComparerManager::Merger::UniquePtr PairwiseCallComparerManager::Merger::create( std::string const& model, appcontext::OptionProcessor const& options ) {
	PairwiseCallComparerManager::Merger::UniquePtr result ;
	if( model == "GenotypeFrequencyTest" ) {
		result.reset(
			new FrequentistTestCallMerger(
				model,
				options.get_value< double >( "-consensus-call-pvalue-threshhold" )
			)
		) ;
	}
	else if( model == "AcceptAll" ) {
		result.reset(
			new AcceptAllCallMerger()
		) ;
	}
	else {
		throw genfile::BadArgumentError( "PairwiseCallComparerManager::Merger::create()", "model=\"" + model + "\"" ) ;
	}
	return result ;
}

