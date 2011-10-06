#ifndef QCTOOL_PAIRWISE_CALL_COMPARER_MANAGER_HPP
#define QCTOOL_PAIRWISE_CALL_COMPARER_MANAGER_HPP

#include <string>
#include <boost/signals2/signal.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "PairwiseCallComparer.hpp"

// This class collects a set of PairwiseCallComparers
// and uses them to compare, pairwise, all the callsets sent to it.
// It sends its output to the given signal.
struct PairwiseCallComparerManager: public genfile::SNPDataSourceProcessor::Callback
{
	virtual ~PairwiseCallComparerManager() {}
	void add_comparer( std::string const& name, PairwiseCallComparer::UniquePtr comparer ) ;
	typedef boost::function< void ( genfile::SNPIdentifyingData const&, std::string const&, std::string const&, std::string const&, std::map< std::string, genfile::VariantEntry > const& ) > ResultCallback ;
	void send_results_to( ResultCallback ) ;
	
	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) ;
	void add_callset_for_snp( genfile::SNPIdentifyingData const& snp, std::string const& name, genfile::SingleSNPGenotypeProbabilities const& probs ) ;
	void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;

private:
	typedef boost::ptr_map< std::string, PairwiseCallComparer > Comparers ;
	Comparers m_comparers ;

	typedef std::map< std::string, genfile::SingleSNPGenotypeProbabilities > Calls ;
	Calls m_calls ;
	genfile::SNPIdentifyingData m_snp ;
	
	typedef boost::signals2::signal< void ( genfile::SNPIdentifyingData const&, std::string const&, std::string const&, std::string const&, std::map< std::string, genfile::VariantEntry > const& ) > ResultSignal ;
	ResultSignal m_result_signal ;
} ;

#endif
