
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_CALLCOMPARERCOMPONENT_CONSENSUS_CALLER_HPP
#define QCTOOL_CALLCOMPARERCOMPONENT_CONSENSUS_CALLER_HPP

#include <boost/signals2/signal.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/VCFFormatSNPDataSink.hpp"
#include "PairwiseCallComparerManager.hpp"

// Given a set of calls passed in by processed_snp
// The consensus caller outputs a 
struct ConsensusCaller: public PairwiseCallComparerManager::MergeClient
{
public:
	typedef std::auto_ptr< ConsensusCaller > UniquePtr ;
	typedef boost::shared_ptr< ConsensusCaller > SharedPtr ;
	static UniquePtr create( std::string const& model ) ;
	static SharedPtr create_shared( std::string const& model ) ;

public:
	ConsensusCaller() ;
public:
	typedef boost::signals2::signal<
		void (
			genfile::SNPIdentifyingData const&,
			Eigen::MatrixXd const&,
			std::map< std::string, std::vector< genfile::VariantEntry > > const&
		)
	> ResultSignal ;

	void send_results_to( ResultSignal::slot_type ) ;
	void send_results(
		genfile::SNPIdentifyingData const& snp,
		Eigen::MatrixXd const& genotypes,
		std::map< std::string, std::vector< genfile::VariantEntry > > const& info
	) ;

	void begin_processing_snps( std::size_t number_of_samples ) ;
	
	void begin_comparisons( genfile::SNPIdentifyingData const& snp ) ;
	virtual void set_result(
		std::string const& comparison,
		std::string const& accepted_calls,
		PairwiseCallComparerManager::Calls const& calls
	) = 0 ;
	void end_comparisons() ;

	std::size_t get_number_of_samples() const { return m_number_of_samples ; }
	genfile::SNPIdentifyingData const& get_snp() const { return m_snp ; }
	std::vector< std::string > const& get_consensus_call_names() const { return m_call_names ; }
private:
	std::vector< std::string > m_call_names ;
	ResultSignal m_result_signal ;
	genfile::SNPIdentifyingData m_snp ;
	std::size_t m_number_of_samples ;
} ;

#endif
