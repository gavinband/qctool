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

void PairwiseCallComparerManager::send_results_to( ResultCallback callback ) {
	m_result_signal.connect( callback ) ;
}


void PairwiseCallComparerManager::set_SNP( genfile::SNPIdentifyingData const& snp ) {
	m_calls.clear() ;
	m_snp = snp ;
}

void PairwiseCallComparerManager::add_calls( std::string const& name, genfile::SingleSNPGenotypeProbabilities const& calls ) {
	m_calls[ name ] = calls ;
}

void PairwiseCallComparerManager::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	m_calls.clear() ;
}

void PairwiseCallComparerManager::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	for( Calls::const_iterator i = m_calls.begin(); i != m_calls.end(); ++i ) {
		Calls::const_iterator j = i ;
		for( ++j; j != m_calls.end(); ++j ) {
			for( Comparers::const_iterator comparer_i = m_comparers.begin(); comparer_i != m_comparers.end(); ++comparer_i ) {
				std::map< std::string, genfile::VariantEntry > const& result = comparer_i->second->compare( i->second, j->second ) ;
				std::map< std::string, genfile::VariantEntry >::const_iterator
					result_i = result.begin(),
					result_end = result.end() ;
				for( ; result_i != result_end; ++result_i ) {
					m_result_signal( snp, i->first, j->first, comparer_i->first, result_i->first, result_i->second ) ;
				}
			}
		}
	}
}

void PairwiseCallComparerManager::end_processing_snps() {
	// nothing to do.
}
