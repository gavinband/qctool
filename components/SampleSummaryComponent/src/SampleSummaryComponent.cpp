
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
#include "components/SampleSummaryComponent/IntensityDistributionComputation.hpp"

void SampleSummaryComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Per-sample computation options" ) ;
    options[ "-sample-stats" ]
		.set_description( "Calculate and output sample-wise statistics." ) ;
    options[ "-sample-stats-file" ]
		.set_description( "Override the default filename of the sample-stats file." )
		.set_takes_single_value() ;
	options.option_implies_option( "-sample-stats", "-s" ) ;
	options.option_implies_option( "-sample-stats-file", "-sample-stats" ) ;
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
	SampleSummaryComputationManager::UniquePtr manager = SampleSummaryComputationManager::create() ;
	std::string filename ;
	if( m_options.check( "-sample-stats-file" ) ) {
		filename = m_options.get_value< std::string >( "-sample-stats-file" ) ;
	}
	else {
		filename = genfile::strip_gen_file_extension_if_present( m_options.get< std::string >( "-g" ) ) + ".qcdb";
	}
	sample_stats::DBOutputter::SharedPtr outputter = sample_stats::DBOutputter::create_shared(
			filename,
			m_options.get< std::string >( "-analysis-name" ),
			m_options.get_values_as_map(),
			m_samples
	) ;
	manager->add_result_callback(
		boost::bind(
			&sample_stats::DBOutputter::operator(),
			outputter,
			_1, _2, _3, _4, _5
		)
	) ;
	manager->add(
		"average intensities",
		"autosomal chromosomes",
		SampleSummaryComputation::UniquePtr( new sample_stats::IntensityDistributionComputation() )
	) ;

	processor.add_callback( genfile::SNPDataSourceProcessor::Callback::UniquePtr( manager.release() ) ) ;
}
