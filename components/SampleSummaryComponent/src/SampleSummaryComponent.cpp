#include <memory>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComponent.hpp"

void SampleSummaryComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Per-sample computation options" ) ;
	options[ "-sample-stats-new" ]
		.set_description( "Compute per-sample statistics." ) ;
}

typedef std::auto_ptr< SampleSummaryComponent > UNiquePtr ;
SampleSummaryComponent::UniquePtr SampleSummaryComponent::create( appcontext::OptionProcessor const& options, genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context ) {
	return SampleSummaryComponent::UniquePtr( new SampleSummaryComponent( options, samples, ui_context )) ;
}

SampleSummaryComponent::SampleSummaryComponent( appcontext::OptionProcessor const& options, genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context ):
	m_options( options ),
	m_samples( samples ),
	m_ui_context( ui_context )
{}

void SampleSummaryComponent::setup( genfile::SNPDataSourceProcessor& processor ) const { 
	// processor.add_callback( ... );
}
