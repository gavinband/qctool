
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_PAIRWISE_CALL_COMPARER_MANAGER_HPP
#define QCTOOL_PAIRWISE_CALL_COMPARER_MANAGER_HPP

#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "components/CallComparerComponent/PairwiseCallComparer.hpp"

// This class collects a set of PairwiseCallComparers
// and uses them to compare, pairwise, all the callsets sent to it.
// It sends its output to the given signal.
struct PairwiseCallComparerManager
{
public:
	typedef std::auto_ptr< PairwiseCallComparerManager > UniquePtr ;
	typedef boost::shared_ptr< PairwiseCallComparerManager > SharedPtr ;
	static UniquePtr create() ;

public:
	PairwiseCallComparerManager() ;
	virtual ~PairwiseCallComparerManager() {}
	void add_comparer( std::string const& name, PairwiseCallComparer::UniquePtr comparer ) ;

	struct Client
	{
		virtual ~Client() {} ;
		virtual void begin_comparisons( genfile::SNPIdentifyingData const& snp ) = 0 ;
		virtual void end_comparisons() = 0 ;
	} ;

	struct ComparisonClient: public virtual Client
	{
		typedef boost::shared_ptr< ComparisonClient > SharedPtr ;
		typedef std::auto_ptr< ComparisonClient > UniquePtr ;

		virtual ~ComparisonClient() {}
		virtual void set_result(
			std::string const& first_callset,
			std::string const& second_callset,
			std::string const& comparison,
			std::string const& comparison_value,
			genfile::VariantEntry const&
		) = 0 ;
	} ;
	
	struct MergeClient: public virtual Client
	{
		typedef boost::shared_ptr< MergeClient > SharedPtr ;
		typedef std::auto_ptr< MergeClient > UniquePtr ;
		virtual ~MergeClient() {}
		virtual void set_result(
			std::string const& comparison,
			std::string const& comparison_value,
			genfile::VariantEntry const&
		) = 0 ;
	} ;
	
	struct Merger: public ComparisonClient {
		typedef std::auto_ptr< Merger > UniquePtr ;
		typedef boost::shared_ptr< Merger > SharedPtr ;
		virtual ~Merger() {}
		virtual std::string get_spec() const = 0 ;
		virtual std::string get_result_as_string() const = 0 ;
	} ;

	void send_comparisons_to( ComparisonClient::SharedPtr ) ;
	void send_merge_to( MergeClient::SharedPtr ) ;
	void set_merger( Merger::UniquePtr ) ;

	void begin_processing_snp( genfile::SNPIdentifyingData const& snp ) ;
	void add_calls( std::string const& name, genfile::SingleSNPGenotypeProbabilities const& calls ) ;
	void end_processing_snp() ;
	
private:
	typedef boost::ptr_map< std::string, PairwiseCallComparer > Comparers ;
	Comparers m_comparers ;
	Merger::UniquePtr m_merger ;
	typedef std::map< std::string, genfile::SingleSNPGenotypeProbabilities > Calls ;
	Calls m_calls ;

	genfile::SNPIdentifyingData m_snp ;

public:
	std::vector< ComparisonClient::SharedPtr > m_comparison_clients ;
	std::vector< MergeClient::SharedPtr > m_merge_clients ;
	void send_comparisons_to_clients(
		std::string const& first_callset,
		std::string const& second_callset,
		std::string const& comparison,
		std::string const& comparison_value,
		genfile::VariantEntry const&
	) ;
	void send_merge_to_clients(
		std::string const& comparison,
		std::string const& variable,
		genfile::VariantEntry const& value
	) ;
} ;


#endif
