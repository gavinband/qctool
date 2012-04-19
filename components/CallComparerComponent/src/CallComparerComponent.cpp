
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <boost/signals2/signal.hpp>
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>


#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "components/CallComparerComponent/PairwiseCallComparer.hpp"
#include "components/CallComparerComponent/PairwiseCallComparerManager.hpp"
#include "components/CallComparerComponent/CallComparerComponent.hpp"
#include "components/CallComparerComponent/FrequentistTestCallMerger.hpp"
#include "components/CallComparerComponent/CallComparerDBOutputter.hpp"
#include "components/CallComparerComponent/CallComparerFileOutputter.hpp"
#include "components/CallComparerComponent/ConsensusCaller.hpp"


CallComparerProcessor::UniquePtr CallComparerProcessor::create( PairwiseCallComparerManager::UniquePtr comparer, std::vector< std::string > const& call_fields ) {
	UniquePtr result(
		new CallComparerProcessor( comparer, call_fields )
	) ;
	return result ;
}

CallComparerProcessor::CallComparerProcessor( PairwiseCallComparerManager::UniquePtr call_comparer, std::vector< std::string > const& call_fields  ):
	m_call_comparer( call_comparer ),
	m_call_fields( call_fields )
{}

void CallComparerProcessor::begin_processing_snps( std::size_t number_of_samples ) {}

void CallComparerProcessor::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	// Add all the calls to the call comparer.
	m_call_comparer->begin_processing_snp( snp ) ;
	genfile::SingleSNPGenotypeProbabilities calls ;
	for( std::size_t i = 0; i < m_call_fields.size(); ++i ) {
		data_reader.get( m_call_fields[i], calls ) ;
		m_call_comparer->add_calls( m_call_fields[i], calls ) ;
	}
	m_call_comparer->end_processing_snp() ;
}

void CallComparerProcessor::end_processing_snps() {}

void CallComparerComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Call comparison options" ) ;
	options[ "-compare-calls" ]
		.set_description( "Compare genotype calls from the given fields of a VCF or other file. "
		 	"The value should be a comma-separated list of genotype call fields in the data.")
		.set_takes_single_value() ;
	options[ "-compare-calls-file" ]
		.set_description( "Name of output file to put call comparisons in." )
		.set_takes_single_value()
		.set_default_value( "call-comparisons.qcdb" ) ;
	options[ "-consensus-call" ]
		.set_description( "Combine calls from several different genotypers into a single consensus callset. "
			"The calls to use are specified using the -compare-calls option. "
			"A consensus call at each SNP is made when no significant difference is found "
			"between N-1 of the N callsets." )
		.set_takes_single_value()
		.set_default_value( "call-comparisons.vcf" ) ;
	options[ "-consensus-strategy" ]
		.set_description( "Strategy to use to combine calls in the consensus set of calls. "
		 	"Currently this must be \"least-missing\"." )
			.set_default_value( "least-missing" ) ;
	options[ "-compare-call-pvalue-threshhold" ]
		.set_description( "Treat calls in call comparison treated as distinct if the p-value of an association test between them "
			"is less than or equal to this value." )
		.set_takes_single_value()
		.set_default_value( 0.001 ) ;

	options.option_implies_option( "-consensus-call", "-s" ) ;
	options.option_implies_option( "-compare-calls", "-s" ) ;
	options.option_implies_option( "-compare-calls", "-consensus-call" ) ;
	options.option_implies_option( "-compare-call-pvalue-threshhold", "-compare-calls" ) ;
}

CallComparerComponent::UniquePtr CallComparerComponent::create(
	genfile::CohortIndividualSource const& samples,
	appcontext::OptionProcessor const& options,
	appcontext::UIContext& ui_context
) {
	CallComparerComponent::UniquePtr result(
		new CallComparerComponent( samples, options, ui_context )
	) ;
	return result ;
}

CallComparerComponent::CallComparerComponent(
	genfile::CohortIndividualSource const& samples,
	appcontext::OptionProcessor const& options,
	appcontext::UIContext& ui_context
):
	m_samples( samples ),
	m_options( options ),
	m_ui_context( ui_context )
{}

void CallComparerComponent::setup( genfile::SNPDataSourceProcessor& processor ) const {
	PairwiseCallComparerManager::UniquePtr manager( PairwiseCallComparerManager::create().release() ) ;

	std::string filename ;
	if( m_options.check( "-compare-calls-file" )) {
		filename = m_options.get< std::string >( "-compare-calls-file" ) ;
	} else {
		filename = genfile::strip_gen_file_extension_if_present( m_options.get< std::string >( "-g" ) ) + ".qcdb";
	}
	if( m_options.check( "-nodb" ) ) {
		CallComparerFileOutputter::SharedPtr outputter = CallComparerFileOutputter::create_shared( filename, m_options.get< std::string >( "-analysis-name" ) ) ;
		manager->send_comparisons_to( outputter ) ;
		manager->send_merge_to( outputter ) ;
	} else {
		CallComparerDBOutputter::SharedPtr outputter = CallComparerDBOutputter::create_shared( filename, m_options.get< std::string >( "-analysis-name" ) ) ;
		manager->send_comparisons_to( outputter ) ;
		manager->send_merge_to( outputter ) ;
	}
	
	ConsensusCaller::SharedPtr consensus_caller(
		new ConsensusCaller(
			genfile::SNPDataSink::create(
				m_options.get< std::string >( "-consensus-call" )
			)
		)
	) ;

	manager->send_merge_to( consensus_caller ) ;
	
	manager->add_comparer(
		"AlleleFrequencyTestCallComparer",
		PairwiseCallComparer::create( "AlleleFrequencyTestCallComparer" )
	) ;

	manager->set_merger(
		PairwiseCallComparerManager::Merger::UniquePtr(
			new FrequentistTestCallMerger(
				"AlleleFrequencyTestCallComparer",
				m_options.get< double >( "-compare-calls-pvalue-threshhold" )
			)
		)
	) ;

	CallComparerProcessor::UniquePtr ccc = CallComparerProcessor::create(
		manager,
		genfile::string_utils::split_and_strip_discarding_empty_entries(
			m_options.get_value< std::string >( "-compare-calls" ),
			",",
			" \t"
		)
	) ;
	
	processor.add_callback( genfile::SNPDataSourceProcessor::Callback::UniquePtr( ccc.release() ) ) ;
	processor.add_callback( *consensus_caller ) ;
}

