
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_HPP

#include <string>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <boost/signals2/signal.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputationManager.hpp"
#include "qcdb/Storage.hpp"

struct SNPSummaryComponent: public boost::noncopyable
{
public:
	typedef std::auto_ptr< SNPSummaryComponent > UniquePtr ;

public:
	static void declare_options( appcontext::OptionProcessor& options ) ;
	static bool is_needed( appcontext::OptionProcessor const& options ) ;

	UniquePtr create(
		genfile::CohortIndividualSource const& samples,
		appcontext::OptionProcessor const& options,
		appcontext::UIContext& ui_context
	) ;
		
	SNPSummaryComponent(
		genfile::CohortIndividualSource const& samples,
		appcontext::OptionProcessor const& options,
		appcontext::UIContext& ui_context
	) ;

	void setup(
		genfile::SNPDataSourceProcessor& processor,
		qcdb::Storage::SharedPtr storage
	) ;
	qcdb::Storage::SharedPtr get_storage() const ;

private:
	stats::SNPSummaryComputationManager::UniquePtr create_manager( qcdb::Storage::SharedPtr ) ;
	void add_computations( stats::SNPSummaryComputationManager& manager, qcdb::Storage::SharedPtr ) const ;
	stats::SNPSummaryComputation::UniquePtr create_computation( std::string const& name ) const ;

private:
	genfile::CohortIndividualSource const& m_samples ;
	appcontext::OptionProcessor const& m_options ;
	appcontext::UIContext& m_ui_context ;
	qcdb::Storage::SharedPtr m_storage ;
} ;

#endif
