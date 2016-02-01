
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
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "components/SNPSummaryComponent/PairwiseCallComparer.hpp"

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
	typedef std::map< std::string, Eigen::MatrixXd > Calls ;

public:
	PairwiseCallComparerManager() ;
	virtual ~PairwiseCallComparerManager() {}
	void add_comparer( std::string const& name, PairwiseCallComparer::UniquePtr comparer ) ;

	struct Client
	{
		virtual ~Client() {} ;
		virtual void begin_comparisons( genfile::VariantIdentifyingData const& snp ) = 0 ;
		virtual void end_comparisons() = 0 ;
	} ;

	struct ComparisonClient: public Client
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
	
	struct MergeClient: public Client
	{
		typedef boost::shared_ptr< MergeClient > SharedPtr ;
		typedef std::auto_ptr< MergeClient > UniquePtr ;
		virtual ~MergeClient() {}
		virtual void begin_processing_snps( std::size_t number_of_samples ) = 0 ;
		virtual void set_result(
			std::string const& comparison_method,
			std::string const& accepted_calls,
			PairwiseCallComparerManager::Calls const&
		) = 0 ;
	} ;
	
	struct Merger: public ComparisonClient {
		typedef std::auto_ptr< Merger > UniquePtr ;
		typedef boost::shared_ptr< Merger > SharedPtr ;
		static UniquePtr create( std::string const&, appcontext::OptionProcessor const& ) ;

		virtual ~Merger() {}
		virtual void add_callset( std::string const& ) = 0 ;

		virtual std::string get_spec() const = 0 ;
		virtual std::string get_result_as_string() const = 0 ;
	} ;



	void send_comparisons_to( ComparisonClient::SharedPtr ) ;
	void send_merge_to( MergeClient::SharedPtr ) ;
	void set_merger( Merger::UniquePtr ) ;

	void begin_processing_snps( std::size_t number_of_samples ) ;
	void begin_processing_snp( genfile::VariantIdentifyingData const& snp ) ;
	void add_calls( std::string const& name, Eigen::MatrixXd const& calls ) ;
	void end_processing_snp() ;
	
private:
	typedef boost::ptr_map< std::string, PairwiseCallComparer > Comparers ;
	Comparers m_comparers ;
	Merger::UniquePtr m_merger ;
	Calls m_calls ;

	genfile::VariantIdentifyingData m_snp ;

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
		std::string const& accepted_calls,
		Calls const& calls
	) ;
} ;


#endif
