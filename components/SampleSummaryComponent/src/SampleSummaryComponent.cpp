
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
#include "components/SampleSummaryComponent/DBOutputter.hpp"
#include "components/SampleSummaryComponent/FlatTableDBOutputter.hpp"
#include "components/SampleSummaryComponent/IntensityDistributionComputation.hpp"
#include "components/SampleSummaryComponent/MissingnessHeterozygosityComputation.hpp"
#include "components/SampleSummaryComponent/RiskScoreComputation.hpp"

void SampleSummaryComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Per-sample computation options" ) ;
    options[ "-sample-stats" ]
		.set_description( "Calculate and output sample-wise statistics." ) ;
    options[ "-sample-stats-table-name" ]
		.set_description( "Set the table name for output of sample statistics when using database output."
		 	" An additional view named <tablename>View will also be created." )
		.set_takes_single_value()
		.set_default_value( "SampleData" )
	;
	options[ "-risk-score" ]
		.set_description( "Compute a risk score for each sample based on a specified file of additive and heterozygote effect sizes." )
		.set_takes_single_value() ;
	
	options.option_implies_option( "-sample-stats", "-o" ) ;
}

SampleSummaryComponent::UniquePtr SampleSummaryComponent::create( appcontext::OptionProcessor const& options, genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context, qcdb::Storage::SharedPtr storage ) {
	return SampleSummaryComponent::UniquePtr( new SampleSummaryComponent( options, samples, ui_context, storage )) ;
}

SampleSummaryComponent::SampleSummaryComponent( appcontext::OptionProcessor const& options, genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context, qcdb::Storage::SharedPtr ):
	m_options( options ),
	m_samples( samples ),
	m_ui_context( ui_context )
{}

void SampleSummaryComponent::setup( genfile::SNPDataSourceProcessor& processor ) const {
	SampleSummaryComputationManager::UniquePtr manager = SampleSummaryComputationManager::create() ;
	std::string filename ;
	if( m_options.check( "-o" ) ) {
		filename = m_options.get_value< std::string >( "-o" ) ;
	}
	else {
		filename = genfile::strip_gen_file_extension_if_present( m_options.get< std::string >( "-g" ) ) + ".qcdb";
	}
	
	{
		sample_stats::FlatTableDBOutputter::SharedPtr outputter = sample_stats::FlatTableDBOutputter::create_shared(
			filename,
			m_options.get< std::string >( "-analysis-name" ),
			m_options.get< std::string >( "-analysis-description" ) + " (sample stats)",
			m_options.get_values_as_map(),
			m_samples
		) ;

		if( m_options.check( "-table-name" ) ) {
			outputter->set_table_name( m_options.get< std::string >( "-table-name" ) + "SampleData" ) ;
		}
	
		manager->send_output_to( sample_stats::SampleStorage::SharedPtr( outputter )) ;
	}
	
	if( m_options.check( "-sample-stats" )) {
		manager->add(
			"average intensities",
			"all chromosomes",
			SampleSummaryComputation::UniquePtr( new sample_stats::IntensityDistributionComputation() )
		) ;

		manager->add(
			"missingness/heterozyosity",
			"all chromosomes",
			SampleSummaryComputation::UniquePtr( new sample_stats::MissingnessHeterozygosityComputation() )
		) ;
	}
	
	if( m_options.check( "-risk-score" )) {
		sample_stats::RiskScoreComputation::UniquePtr computation(
			new sample_stats::RiskScoreComputation(
				m_samples,
				genfile::SNPIdentifyingData::CompareFields( m_options.get_value< std::string >( "-snp-match-fields" ) )
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

	m_ui_context.logger() << "SampleSummaryComponent: the following components are in place:\n" << manager->get_summary( "  " ) << "\n" ;

	processor.add_callback( genfile::SNPDataSourceProcessor::Callback::UniquePtr( manager.release() ) ) ;
}
