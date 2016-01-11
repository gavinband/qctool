
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <memory>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComponent.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComputationManager.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"
#include "components/SampleSummaryComponent/SampleStorage.hpp"
#include "components/SampleSummaryComponent/IntensityDistributionComputation.hpp"
#include "components/SampleSummaryComponent/MissingnessHeterozygosityComputation.hpp"
#include "components/SampleSummaryComponent/RiskScoreComputation.hpp"

void SampleSummaryComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Per-sample computation options" ) ;
    options[ "-sample-stats" ]
		.set_description( "Calculate and output sample-wise statistics." ) ;
	options[ "-risk-score" ]
		.set_description( "Compute a risk score for each sample based on a specified file of additive and heterozygote effect sizes." )
		.set_takes_single_value() ;
	
	options.option_implies_option( "-sample-stats", "-osample" ) ;
}

bool SampleSummaryComponent::is_needed( appcontext::OptionProcessor const& options ) {
	return options.check( "-sample-stats" )
		|| options.check( "-risk-score" )
	;
}

SampleSummaryComponent::UniquePtr SampleSummaryComponent::create(
	appcontext::OptionProcessor const& options,
	genfile::CohortIndividualSource const& samples,
	appcontext::UIContext& ui_context
) {
	return SampleSummaryComponent::UniquePtr( new SampleSummaryComponent( options, samples, ui_context )) ;
}

SampleSummaryComponent::SampleSummaryComponent(
	appcontext::OptionProcessor const& options,
	genfile::CohortIndividualSource const& samples,
	appcontext::UIContext& ui_context
):
	m_options( options ),
	m_samples( samples ),
	m_ui_context( ui_context )
{}

void SampleSummaryComponent::setup(
	genfile::SNPDataSourceProcessor& processor,
	sample_stats::SampleStorage::SharedPtr storage
) const {
	SampleSummaryComputationManager::UniquePtr manager = SampleSummaryComputationManager::create() ;
	
	if( m_options.check( "-sample-stats" )) {
		manager->add(
			"missingness/heterozyosity",
			"all chromosomes",
			SampleSummaryComputation::UniquePtr( new sample_stats::MissingnessHeterozygosityComputation() )
		) ;

		manager->add(
			"average intensities",
			"all chromosomes",
			SampleSummaryComputation::UniquePtr( new sample_stats::IntensityDistributionComputation() )
		) ;
	}
	
	if( m_options.check( "-risk-score" )) {
		sample_stats::RiskScoreComputation::UniquePtr computation(
			new sample_stats::RiskScoreComputation(
				m_samples,
				genfile::VariantIdentifyingData::CompareFields( m_options.get_value< std::string >( "-snp-match-fields" ) )
			)
		) ;

		{
			appcontext::UIContext::ProgressContext progress_context = m_ui_context.get_progress_context( "Loading effects from \"" + m_options.get_value< std::string >( "-risk-score" ) + "\"" ) ;
			statfile::BuiltInTypeStatSource::UniquePtr source = statfile::BuiltInTypeStatSource::open( genfile::wildcard::find_files_by_chromosome( m_options.get_value< std::string >( "-risk-score" ) ) ) ;
			computation->add_effects( *source, progress_context ) ;
		}

		manager->add(
			"risk_score",
			"all chromosomes",
			SampleSummaryComputation::UniquePtr(
				computation.release()
			)
		) ;
	}

	manager->send_output_to( sample_stats::SampleStorage::SharedPtr( storage )) ;

	m_ui_context.logger() << "SampleSummaryComponent: the following components are in place:\n" << manager->get_summary( "  " ) << "\n" ;
	processor.add_callback( genfile::SNPDataSourceProcessor::Callback::UniquePtr( manager.release() ) ) ;
}
