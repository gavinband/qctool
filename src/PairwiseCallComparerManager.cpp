#include <string>
#include <boost/signals2/signal.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "PairwiseCallComparer.hpp"
#include "PairwiseCallComparerManager.hpp"

void PairwiseCallComparerManager::add_comparer( std::string const& name, PairwiseCallComparer::UniquePtr comparer ) {
	m_comparers.insert( name, comparer ) ;
}

void PairwiseCallComparerManager::send_results_to( ResultCallback callback ) {
	m_result_signal.connect( callback ) ;
}

void PairwiseCallComparerManager::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	// do nothing.
	m_calls.clear() ;
}

void PairwiseCallComparerManager::add_callset_for_snp(
	genfile::SNPIdentifyingData const& snp,
	std::string const& name,
	genfile::SingleSNPGenotypeProbabilities const& probs
) {
	m_calls[ name ] = probs ;
}

void PairwiseCallComparerManager::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	for( Calls::const_iterator i = m_calls.begin(); i != m_calls.end(); ++i ) {
		Calls::const_iterator j = i ;
		for( ++j; j != m_calls.end(); ++j ) {
			for( Comparers::const_iterator comparer_i = m_comparers.begin(); comparer_i != m_comparers.end(); ++comparer_i ) {
				m_result_signal( snp, i->first, j->first, comparer_i->first, comparer_i->second->compare( i->second, j->second )) ;
			}
		}
	}
	m_calls.clear() ;
}

void PairwiseCallComparerManager::end_processing_snps() {
	// nothing to do.
}
