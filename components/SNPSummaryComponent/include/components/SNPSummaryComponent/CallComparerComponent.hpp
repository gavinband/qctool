
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
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputationManager.hpp"
#include "qcdb/Storage.hpp"

struct CallComparerProcessor: public SNPSummaryComputation
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

		void operator()( SNPIdentifyingData const&, Genotypes const&, SampleSexes const&, genfile::VariantDataReader&, ResultCallback ) ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;

	private:
		PairwiseCallComparerManager::UniquePtr m_call_comparer ;
		std::vector< std::string > m_call_fields ;
		bool m_begun ;
		PairwiseCallComparerManager::MergeClient::UniquePtr m_consensus_caller ;
		Eigen::MatrixXd m_calls ;
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
	
	void setup( SNPSummaryComputationManager&, qcdb::Storage::SharedPtr ) const ;

private:
	genfile::CohortIndividualSource const& m_samples ;
	appcontext::OptionProcessor const& m_options ;
	appcontext::UIContext& m_ui_context ; 
} ;


#endif
