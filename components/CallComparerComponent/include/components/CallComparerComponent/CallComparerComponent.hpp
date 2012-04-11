
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_CALL_COMPARER_COMPONENT_HPP
#define QCTOOL_CALL_COMPARER_COMPONENT_HPP

#include <boost/noncopyable.hpp>
#include "genfile/CohortIndividualSource.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "PairwiseCallComparerManager.hpp"

struct CallComparerProcessor: public genfile::SNPDataSourceProcessor::Callback, public virtual boost::noncopyable
{
public:
	typedef std::auto_ptr< CallComparerProcessor > UniquePtr ;
	static UniquePtr create(
		PairwiseCallComparerManager::UniquePtr call_comparer,
		std::vector< std::string > const& call_fields
	) ;

public:
	CallComparerProcessor(
		PairwiseCallComparerManager::UniquePtr call_comparer,
		std::vector< std::string > const& call_fields
	) ;

	public:
		void begin_processing_snps( std::size_t number_of_samples ) ;
		void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& data_reader ) ;
		void end_processing_snps() ;

		PairwiseCallComparerManager& get_comparer() { return *m_call_comparer ; }
		
	private:
		PairwiseCallComparerManager::UniquePtr m_call_comparer ;
		std::vector< std::string > m_call_fields ;
		PairwiseCallComparerManager::MergeClient::UniquePtr m_consensus_caller ;
} ;

struct CallComparerComponent: public boost::noncopyable
{
	static void declare_options( appcontext::OptionProcessor& options ) ;

	typedef std::auto_ptr< CallComparerComponent > UniquePtr ;
	static UniquePtr create(
		genfile::CohortIndividualSource const& samples,
		appcontext::OptionProcessor const& options,
		appcontext::UIContext& ui_context
	) ;

public:
	CallComparerComponent(
		genfile::CohortIndividualSource const& m_samples,
		appcontext::OptionProcessor const& options,
		appcontext::UIContext& ui_context
	) ;
	
	void setup( genfile::SNPDataSourceProcessor& processor ) const ;

private:
	genfile::CohortIndividualSource const& m_samples ;
	appcontext::OptionProcessor const& m_options ;
	appcontext::UIContext& m_ui_context ;
} ;


#endif
