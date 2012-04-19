
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComponent.hpp"
#include "../../qctool_version_autogenerated.hpp"
#include "components/SNPSummaryComponent/DBOutputter.hpp"
#include "components/SNPSummaryComponent/FileOutputter.hpp"
#include "components/SNPSummaryComponent/AssociationTest.hpp"
#include "components/SNPSummaryComponent/AncestralAlleleAnnotation.hpp"

void SNPSummaryComputationManager::add_computation( std::string const& name, SNPSummaryComputation::UniquePtr computation ) {
	m_computations.insert( name, computation ) ;
}

void SNPSummaryComputationManager::add_result_callback( ResultCallback callback ) {
	m_result_signal.connect( callback ) ;
}

void SNPSummaryComputationManager::begin_processing_snps( std::size_t number_of_samples ) {
	m_snp_index = 0 ;
	m_genotypes.resize( number_of_samples, 3 ) ;
}

void SNPSummaryComputationManager::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	{
		genfile::vcf::GenotypeSetter< Eigen::MatrixBase< SNPSummaryComputation::Genotypes > > setter( m_genotypes ) ;
		data_reader.get( "genotypes", setter ) ;
	}
	Computations::iterator i = m_computations.begin(), end_i = m_computations.end() ;
	for( ; i != end_i; ++i ) {
		i->second->operator()(
			snp,
			m_genotypes,
			data_reader,
			boost::bind(
				boost::ref( m_result_signal ),
				m_snp_index,
				snp,
				i->first,
				_1,
				_2
			)
		) ;
	}
	++m_snp_index ;
}

void SNPSummaryComputationManager::end_processing_snps() {}

std::string SNPSummaryComputationManager::get_summary( std::string const& prefix, std::size_t column_width ) {
	std::string result ;
	Computations::const_iterator i = m_computations.begin(), end_i = m_computations.end() ;
	for( ; i != end_i; ++i ) {
		result += i->second->get_summary( prefix, column_width ) + "\n";
	}
	return result ;
}

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
	
	options[ "-annotate" ]
		.set_description( "Specify a FASTA-formatted file containing ancestral alleles to annotate variants with." )
		.set_takes_single_value() ;
	
	options.option_implies_option( "-snp-stats", "-g" ) ;
	options.option_implies_option( "-annotate", "-g" ) ;
	options.option_implies_option( "-test", "-g" ) ;
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
	SNPSummaryComputationManager::UniquePtr manager( new SNPSummaryComputationManager() ) ;
	using genfile::string_utils::to_string ;
	
	std::string filename ;
	if( m_options.check( "-snp-stats" )) {
		filename = m_options.get_value< std::string >( "-snp-stats" ) ;
	}
	else {
		std::vector< std::string > filenames = m_options.get_values< std::string >( "-g" ) ;
		if( filenames.size() == 1 ) {
			filename = genfile::strip_gen_file_extension_if_present( filenames[0] ) + ".snp-stats.txt" ;
		} else {
			filename = "qctool_cohort_1-" + to_string( filenames.size() ) + ".snp-stats.txt" ;
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
					phenotype,
					covariates,
					m_samples,
					m_options
				)
			) ;
		}
	}
	if( m_options.check( "-annotate" )) {
		appcontext::UIContext::ProgressContext progress = m_ui_context.get_progress_context( "Loading FASTA annotation" ) ;
		manager.add_computation(
			"ancestral_alleles",
			SNPSummaryComputation::UniquePtr(
				new AncestralAlleleAnnotation( m_options.get< std::string >( "-annotate" ), progress )
			)
		) ;
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
