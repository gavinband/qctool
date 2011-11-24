#include <boost/math/distributions/chi_squared.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/VariantEntry.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "snptest/PerSnpFrequentistTest.hpp"
#include "integration/Derivative.hpp"
#include "integration/NewtonRaphson.hpp"
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
	options[ "-mimic-snptest" ]
		.set_description(
			"Configure things so the tests agree as far as possible with those of SNPTEST (v2.2.0)."
			" This affects how qctool deals with NULL genotype calls."
		) ;

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
	m_phenotypes( get_phenotypes( samples, m_options.get_value< std::string >( "-test" ))),
	m_indices_of_samples_to_include( get_indices_of_samples_to_include( m_phenotypes, m_samples )),
	m_phenotype_values( get_phenotype_values( m_phenotypes, m_samples, m_indices_of_samples_to_include ))
{
	assert( m_phenotypes.size() > 0 ) ;
}

std::vector< std::string > AssociationTester::get_phenotypes(
	genfile::CohortIndividualSource const& samples,
	std::string const& phenotype_spec
) const {
	std::vector< std::string > phenotypes = genfile::string_utils::split_and_strip( phenotype_spec, ",", " \t\n" ) ;
	genfile::CohortIndividualSource::ColumnSpec column_spec = samples.get_column_spec() ;
	for( std::size_t i = 0; i < phenotypes.size(); ++i ) {
		column_spec.find_column( phenotypes[i] ) ;
	}
	return phenotypes ;
}

std::vector< std::vector< std::size_t > > AssociationTester::get_indices_of_samples_to_include(
	std::vector< std::string > const& phenotypes,
	genfile::CohortIndividualSource const& samples
) const {
	std::vector< std::vector< std::size_t > > result( phenotypes.size() ) ;
	for( std::size_t i = 0; i < phenotypes.size(); ++i ) {
		for( std::size_t j = 0; j < samples.get_number_of_individuals(); ++j ) {
			if( !samples.get_entry( j, phenotypes[i] ).is_missing() ) {
				result[i].push_back( j ) ;
			}
		}
	}
	return result ;
}

std::vector< AssociationTester::Vector > AssociationTester::get_phenotype_values(
	std::vector< std::string > const& phenotypes,
	genfile::CohortIndividualSource const& samples,
	std::vector< std::vector< std::size_t > > indices_of_samples_to_include
) const {
	assert( indices_of_samples_to_include.size() == phenotypes.size() ) ;
	std::vector< AssociationTester::Vector > result( phenotypes.size() ) ;
	for( std::size_t i = 0; i < phenotypes.size(); ++i ) {
		result[i] = Vector::Zero( indices_of_samples_to_include[i].size() ) ;
		for( std::size_t j = 0; j < indices_of_samples_to_include[i].size(); ++j ) {
			result[i]( j ) = double( samples.get_entry( indices_of_samples_to_include[i][j], phenotypes[i] ).as< int >() ) ;
		}
	}
	return result ;
}

void AssociationTester::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	assert( number_of_samples == m_samples.get_number_of_individuals() ) ;
	m_sink = statfile::BuiltInTypeStatSink::open( m_options.get_value< std::string >( "-ot" )) ;
	m_sink->write_metadata( "This file written by qctool, " + appcontext::get_current_time_as_string() + ".\n" ) ;
	(*m_sink) | "SNPID" | "rsid" | "chromosome" | "position" | "alleleA" | "alleleB" ;
	for( std::size_t i = 0; i < m_phenotypes.size(); ++i ) {
		(*m_sink) | ( m_phenotypes[i] + ".p-value" ) | ( m_phenotypes[i] + ".beta" ) | ( m_phenotypes[i] + ".standard_error" ) ;
		(*m_sink) | ( m_phenotypes[i] + ".sigma_11" ) | ( m_phenotypes[i] + ".sigma_12" ) | ( m_phenotypes[i] + ".sigma_21" ) | ( m_phenotypes[i] + ".sigma_22" ) ;
	}
}

void AssociationTester::processed_snp(
	SNPIdentifyingData const& id_data,
	genfile::VariantDataReader& data_reader
) {
	SingleSNPGenotypeProbabilities genotypes( m_samples.get_number_of_individuals() ) ;
	data_reader.get( "genotypes", genotypes ) ;
	processed_snp(
		id_data,
		genotypes
		) ;
}

void AssociationTester::processed_snp(
	SNPIdentifyingData const& id_data,
	SingleSNPGenotypeProbabilities const& raw_genotypes
) {
	assert( m_sink.get() ) ;
	(*m_sink)
		<< id_data.get_SNPID()
		<< id_data.get_rsid()
		<< std::string( id_data.get_position().chromosome() )
		<< id_data.get_position().position()
		<< id_data.get_first_allele()
		<< id_data.get_second_allele()
	;
	
	for( std::size_t i = 0; i < m_phenotypes.size(); ++i ) {
		snptest::FinitelySupportedFunctionSet genotypes = get_genotype_matrix( raw_genotypes, m_indices_of_samples_to_include[i] ) ;
		
		snptest::PerSnpFrequentistTest::UniquePtr test = snptest::PerSnpFrequentistTest::create(
			id_data,
			m_options
		) ;

		snptest::PerSnpFrequentistTest::Results results = test->test(
			m_phenotype_values[i],
			m_covariate_values,
			id_data,
			genotypes
		) ;

		(*m_sink)
			<< results.p_value
			<< results.beta
			<< results.standard_error
			<< results.variance_covariance(0,0)
			<< results.variance_covariance(0,1)
			<< results.variance_covariance(1,0)
			<< results.variance_covariance(1,1)
		;
	}
	(*m_sink) << statfile::end_row() ;
}

snptest::FinitelySupportedFunctionSet AssociationTester::get_genotype_matrix(
	genfile::SingleSNPGenotypeProbabilities const& genotypes,
	std::vector< std::size_t > const& indices_of_samples_to_include
 ) const {
	std::size_t const N = indices_of_samples_to_include.size() ;
	double const NaN = std::numeric_limits< double >::quiet_NaN() ;
	snptest::FinitelySupportedFunctionSet result(
		Vector::Constant( 3, NaN ),
		Matrix::Constant( N, 3, NaN )
	) ;

	for( std::size_t i = 0; i < N; ++i ) {
		for( std::size_t g = 0; g < 3; ++g ) {
			result.get_values( i )( g ) = genotypes( indices_of_samples_to_include[ i ], g ) ;
		}
	}
	// mean-centre the genotypes.
	double mean_genotype = ( result.get_values().col(1) + 2.0 * result.get_values().col(2) ).sum() / N ;
	for( std::size_t g = 0; g < 3; ++g ) {
		result.get_support( g ) = double( g ) - mean_genotype ;
	}
	
	return result ;
}

void AssociationTester::end_processing_snps() {
	m_sink.reset() ;
	m_ui_context.logger() << "AssociationTester: finished processing SNPs.\n" ;
}

