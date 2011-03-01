#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/string_utils.hpp"
#include "AssociationTester.hpp"

void AssociationTester::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Association test options" ) ;
	options[ "-test" ]
		.set_description( "Perform an association test on the given phenotype." )
		.set_takes_single_value() ;
	options[ "-covariates" ]
		.set_description( "Specify a comma-separated list of covariates to use in the association test." )
		.set_takes_single_value()
		.set_default_value( "" ) ;
	options[ "-ot" ]
		.set_description( "Override the default name of the association test output file." )
		.set_takes_single_value()
		.set_default_value( "qctool.test" ) ;

	options.option_implies_option( "-ot", "-test" ) ;
	options.option_implies_option( "-covariates", "-test" ) ;
	options.option_implies_option( "-test", "-s" ) ;
}

AssociationTester::AssociationTester(
	appcontext::OptionProcessor const& options,
	genfile::CohortIndividualSource const& samples,
	appcontext::UIContext& ui_context
):
	m_options( options ),
	m_samples( samples ),
	m_ui_context( ui_context ),
	m_phenotypes( genfile::string_utils::split_and_strip( m_options.get_value< std::string >( "-test" ), ",", " \n\t" ))
{
	assert( m_phenotypes.size() > 0 ) ;
	for( std::size_t i = 0; i < m_phenotypes.size(); ++i ) {
		assert( m_samples.check_for_column( m_phenotypes[i] )) ;
	}
}

void AssociationTester::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	assert( number_of_samples == m_samples.get_number_of_individuals() ) ;
	m_sink = statfile::BuiltInTypeStatSink::open( m_options.get_value< std::string >( "-ot" )) ;
	(*m_sink) | "SNPID" | "rsid" | "chromosome" | "position" | "alleleA" | "alleleB" | "p-value" | "beta" | "standard.error" ;
}

void AssociationTester::processed_snp( SNPIdentifyingData const& id_data, SingleSNPGenotypeProbabilities const& genotypes ) {
	assert( m_sink.get() ) ;
	(*m_sink)
		<< id_data.get_SNPID()
		<< id_data.get_rsid()
		<< std::string( id_data.get_position().chromosome() )
		<< id_data.get_position().position()
		<< id_data.get_first_allele()
		<< id_data.get_second_allele()
		<< "NA"
		<< "NA"
		<< "NA"
		<< statfile::end_row() ;
}

void AssociationTester::end_processing_snps() {
	m_sink.reset() ;
	m_ui_context.logger() << "AssociationTester: finished processing SNPs.\n" ;
}

