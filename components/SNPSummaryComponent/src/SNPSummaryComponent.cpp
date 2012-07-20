
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/Error.hpp"
#include "../../qctool_version_autogenerated.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComponent.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputationManager.hpp"
#include "components/SNPSummaryComponent/DBOutputter.hpp"
#include "components/SNPSummaryComponent/FileOutputter.hpp"
#include "components/SNPSummaryComponent/AssociationTest.hpp"
#include "components/SNPSummaryComponent/SequenceAnnotation.hpp"
#include "components/SNPSummaryComponent/DifferentialMissingnessComputation.hpp"
#include "components/SNPSummaryComponent/StratifyingSNPSummaryComputation.hpp"

void SNPSummaryComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "SNP computation options" ) ;
	options[ "-snp-stats" ]
		.set_description( "Calculate and output per-SNP statistics.  This implies that no SNP filtering options are used." ) ;

	options[ "-snp-stats-columns" ]
        .set_description( "Comma-seperated list of extra columns to output in the snp-wise statistics file." )
		.set_takes_single_value()
		.set_default_value( "alleles,HWE,missingness,information" ) ;

	options.declare_group( "Association test options" ) ;
	options[ "-test" ]
		.set_description( "Perform an association test on the given phenotype." )
		.set_takes_single_value() ;
	options[ "-covariates" ]
		.set_description( "Specify a comma-separated list of covariates to use in the association test." )
		.set_takes_single_value()
		.set_default_value( "" ) ;
	options[ "-no-X-inactivation" ]
		.set_description( "Specify that X chromosome inactivation in females should not be modelled in the association test. "
			"If this option is specified, females have twice the maximum exposure that males do." )
		.set_takes_single_value()
		.set_default_value( "" ) ;
	
	options[ "-stratify" ]
		.set_description( "Compute all SNP summary statistics seperately for each level of the given variable in the sample file." )
		.set_takes_single_value() ;

	options[ "-differential" ]
		.set_description( "Test for differences in SNP summary statistics between the categories of the given variable.  Currently a test for differential missingness is performed." )
		.set_takes_single_value() ;
	
	options.declare_group( "Sequence annotation options" ) ;
	options[ "-annotate-ancestral" ]
		.set_description( "Specify a FASTA-formatted file containing ancestral alleles to annotate variants with." )
		.set_takes_single_value() ;
	options[ "-annotate-reference" ]
		.set_description( "Specify a FASTA-formatted file containing ancestral alleles to annotate variants with." )
		.set_takes_single_value() ;
	options[ "-flanking" ]
		.set_description( "Specify that flanking sequence annotations [ pos - a, pos + b ] should be output when using "
			"-annotate-reference and -annotate-ancestral" )
		.set_takes_values( 2 )
		.set_minimum_multiplicity( 0 )
		.set_maximum_multiplicity( 1 ) ;

	
	options.option_implies_option( "-snp-stats", "-g" ) ;
	options.option_implies_option( "-annotate", "-g" ) ;
	options.option_implies_option( "-flanking", "-annotate" ) ;
	options.option_implies_option( "-test", "-g" ) ;
	options.option_implies_option( "-test", "-s" ) ;
	options.option_implies_option( "-stratify", "-s" ) ;
}

SNPSummaryComponent::SNPSummaryComponent(
	genfile::CohortIndividualSource const& samples,
	appcontext::OptionProcessor const& options,
	appcontext::UIContext& ui_context
):
	m_samples( samples ),
	m_options( options ),
	m_ui_context( ui_context )
{}

genfile::SNPDataSourceProcessor::Callback::UniquePtr SNPSummaryComponent::create() const {
	genfile::SNPDataSourceProcessor::Callback::UniquePtr result( create_manager().release() ) ;
	return result ;
}

SNPSummaryComputationManager::UniquePtr SNPSummaryComponent::create_manager() const {
	SNPSummaryComputationManager::UniquePtr manager( new SNPSummaryComputationManager( m_samples ) ) ;
	using genfile::string_utils::to_string ;
	
	std::string filename ;
	if( m_options.check( "-odb" )) {
		filename = m_options.get_value< std::string >( "-odb" ) ;
	}
	else {
		std::vector< std::string > filenames = m_options.get_values< std::string >( "-g" ) ;
		if( filenames.size() == 1 ) {
			filename = genfile::strip_gen_file_extension_if_present( filenames[0] ) + ( m_options.check( "-nodb" ) ? ".snp-stats.tsv" : ".qcdb" ) ;
		} else {
			filename = "qctool_cohort_1-" + to_string( filenames.size() ) + ( m_options.check( "-nodb" ) ? ".snp-stats.tsv" : ".qcdb" ) ;
		}
	}

	if( m_options.check( "-nodb" )) {
		manager->add_result_callback(
			boost::bind(
				&FileOutputter::operator(),
				FileOutputter::create_shared( filename ),
				_1, _2, _3, _4, _5
			)
		) ;
	}
	else {
		snp_summary_component::DBOutputter::SharedPtr outputter = snp_summary_component::DBOutputter::create_shared(
			filename,
			m_options.get< std::string >( "-analysis-name" ),
			m_options.get_values_as_map()
		) ;
	
		manager->add_result_callback(
			boost::bind(
				&snp_summary_component::DBOutputter::operator(),
				outputter,
				_1, _2, _3, _4, _5
			)
		) ;
	}
	
	add_computations( *manager ) ;
	return manager ;
}

namespace impl {
	
	typedef std::map< genfile::VariantEntry, std::vector< int > > StrataMembers ;

	StrataMembers compute_strata( genfile::CohortIndividualSource const& samples, std::string const& variable ) {
		StrataMembers result ;
		genfile::CohortIndividualSource::ColumnSpec const spec = samples.get_column_spec() ;
		if( !spec[ variable ].is_discrete() ) {
			throw genfile::BadArgumentError( "void impl::compute_strata()", "variable=\"" + variable + "\"" ) ;
		}

		for( std::size_t i = 0; i < samples.get_number_of_individuals(); ++i ) {
			genfile::VariantEntry const& entry = samples.get_entry( i, variable ) ;
			if( !entry.is_missing() ) {
				result[ samples.get_entry( i, variable ) ].push_back( int( i ) ) ;
			}
		}

		return result ;
	}
}

void SNPSummaryComponent::add_computations( SNPSummaryComputationManager& manager ) const {
	using genfile::string_utils::split_and_strip_discarding_empty_entries ;

	if( m_options.check( "-snp-stats" )) {
		std::vector< std::string > elts = split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "-snp-stats-columns" ), ",", " \t" ) ;
		foreach( std::string const& elt, elts ) {
			manager.add_computation( elt, SNPSummaryComputation::create( elt )) ;
		}
	}

	if( m_options.check( "-test" )) {
		std::vector< std::string > phenotypes = split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "-test" ), ",", " \t" ) ;
		std::vector< std::string > covariates ;
		if( m_options.check( "-covariates" ))  {
			covariates = split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "-covariates" ), ",", " \t" ) ;
		}
		foreach( std::string const& phenotype, phenotypes ) {
			manager.add_computation(
				"association_test",
				AssociationTest::create(
					"autosomal",
					phenotype,
					covariates,
					m_samples,
					m_options
				)
			) ;
			manager.add_computation(
				"X_chromosome_association_test",
				AssociationTest::create(
					"X chromosome",
					phenotype,
					covariates,
					m_samples,
					m_options
				)
			) ;
		}
	}

	if( m_options.check( "-annotate-ancestral" )) {
		appcontext::UIContext::ProgressContext progress = m_ui_context.get_progress_context( "Loading ancestral sequence" ) ;
		SequenceAnnotation::UniquePtr computation(
			new SequenceAnnotation( "ancestral", m_options.get< std::string >( "-annotate-ancestral" ), progress )
		) ;
		
		if( m_options.check( "-flanking" )) {
			std::vector< std::size_t > data = m_options.get_values< std::size_t >( "-flanking" ) ;
			assert( data.size() == 2 ) ;
			computation->set_flanking( data[0], data[1] ) ;
		}

		manager.add_computation(
			"ancestral_sequence",
			SNPSummaryComputation::UniquePtr(
				computation.release()
			)
		) ;
	}

	if( m_options.check( "-annotate-reference" )) {
		appcontext::UIContext::ProgressContext progress = m_ui_context.get_progress_context( "Loading reference sequence" ) ;
		SequenceAnnotation::UniquePtr computation(
			new SequenceAnnotation( "reference", m_options.get< std::string >( "-annotate-reference" ), progress )
		) ;
		
		if( m_options.check( "-flanking" )) {
			std::vector< std::size_t > data = m_options.get_values< std::size_t >( "-flanking" ) ;
			assert( data.size() == 2 ) ;
			computation->set_flanking( data[0], data[1] ) ;
		}
		
		manager.add_computation(
			"reference_sequence",
			SNPSummaryComputation::UniquePtr(
				computation.release()
			)
		) ;
	}

	if( m_options.check( "-differential" )) {
		std::string const variable = m_options.get< std::string >( "-differential" ) ;
		impl::StrataMembers strata = impl::compute_strata( m_samples, variable ) ;
		manager.add_computation(
			"differential_missingness",
			DifferentialMissingnessComputation::create( variable, strata )
		) ;
	}

	if( m_options.check( "-stratify" )) {
		std::string const variable = m_options.get< std::string >( "-stratify" ) ;
		impl::StrataMembers strata = impl::compute_strata( m_samples, variable ) ;
		manager.stratify_by( strata, variable ) ;
	}
	
	m_ui_context.logger() << "SNPSummaryComponent: the following components are in place:\n" << manager.get_summary( "  " ) << "\n" ;
}

SNPSummaryComputation::UniquePtr SNPSummaryComponent::create_computation( std::string const& name ) const {
	if( name != "association_test" ) {
		return SNPSummaryComputation::UniquePtr( SNPSummaryComputation::create( name )) ;
	} else {
		assert(0) ;
	}
}
