#ifndef QCTOOL_CALLCOMPARERCOMPONENT_CONSENSUS_CALLER_HPP
#define QCTOOL_CALLCOMPARERCOMPONENT_CONSENSUS_CALLER_HPP

#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/VCFFormatSNPDataSink.hpp"
#include "PairwiseCallComparerManager.hpp"

// Given a set of calls passed in by processed_snp
// The consensus caller outputs a 
struct ConsensusCaller: public PairwiseCallComparerManager::MergeClient, public genfile::SNPDataSourceProcessor::Callback
{
public:
	typedef std::auto_ptr< ConsensusCaller > UniquePtr ;
	typedef boost::shared_ptr< ConsensusCaller > SharedPtr ;
	ConsensusCaller( genfile::SNPDataSink::UniquePtr sink ) ;
public:
	void begin_processing_snps( std::size_t number_of_samples ) ;
	void begin_comparisons( genfile::SNPIdentifyingData const& snp ) ;
	void set_result(
		std::string const& comparison,
		std::string const& comparison_value,
		genfile::VariantEntry const&
	) ;
	void end_comparisons() ;
	void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() {}
private:
	double const m_call_threshhold ;
	std::vector< std::string > m_call_names ;
	std::vector< Eigen::MatrixXd > m_genotypes ;
	genfile::SNPDataSink::UniquePtr m_sink ;
} ;

#endif
