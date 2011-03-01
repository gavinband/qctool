#include <boost/math/distributions/chi_squared.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/string_utils.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "snptest/SNPTEST2NullModel.hpp"
#include "snptest/SNPTEST2AlternativeModel.hpp"
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
	m_phenotype_values( get_phenotype_values( m_phenotypes, m_samples, m_indices_of_samples_to_include )),
	m_fill_null_genotypes( false ),
	m_chi_squared( 1.0 )
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

void AssociationTester::processed_snp( SNPIdentifyingData const& id_data, SingleSNPGenotypeProbabilities const& genotypes ) {
	assert( m_sink.get() ) ;
	(*m_sink)
		<< id_data.get_SNPID()
		<< id_data.get_rsid()
		<< std::string( id_data.get_position().chromosome() )
		<< id_data.get_position().position()
		<< id_data.get_first_allele()
		<< id_data.get_second_allele() ;
	for( std::size_t i = 0; i < m_phenotypes.size(); ++i ) {
		FrequentistTestResults results = get_frequentist_test_results(
			m_phenotype_values[i],
			genotypes,
			m_indices_of_samples_to_include[i]
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

AssociationTester::FrequentistTestResults AssociationTester::get_frequentist_test_results(
	Vector const& phenotype_values,
	SingleSNPGenotypeProbabilities const& all_genotypes,
	std::vector< std::size_t > const& indices_of_samples_to_include
) const {
	assert( std::size_t( phenotype_values.size() ) == indices_of_samples_to_include.size() ) ;
	Matrix genotypes( phenotype_values.size(), 3 ) ;
	std::vector< double > genotype_levels(3) ;
	
	{
		for( int i = 0; i < phenotype_values.size(); ++i ) {
			for( std::size_t g = 0; g < 3; ++g ) {
				genotypes( i, g ) = all_genotypes( indices_of_samples_to_include[i], g ) ;
			}
		}
		// We now deal with missing genotypes (= null calls).
		// We fill in any null calls from the estimated allele frequency.
		// The allele frequency estimate is a maximum likelihood estimate,
		// but for the likelihood which treates the calls as the actual data.
		// In practice I guess this is ok.
		double I = genotypes.sum() ;
		double theta = ( genotypes.col(1) + 2.0 * genotypes.col(2) ).sum() / (2.0 * I ) ;
		std::size_t const N = genotypes.rows() ;

		Vector expected_genotypes( 3 ) ;
		expected_genotypes
			<< ( 1.0 - theta ) * ( 1.0 - theta ),
			2.0 * theta * ( 1.0 - theta ),
			theta * theta ;
		if( m_fill_null_genotypes ) {
			for( std::size_t i = 0; i < N; ++i ) {
				double sum = genotypes.row(i).sum() ;
				if( sum > 1.0 ) {
					// assume this is due to rounding error.
					genotypes.row(i) /= sum ;
				}
				else {
					// fill in with expected genotype.
					genotypes.row(i) += ( 1.0 - sum ) * expected_genotypes ;
				}
			}
		}
		double mean_genotype = ( genotypes.col(1) + 2.0 * genotypes.col(2) ).sum() / N ;
		for( std::size_t g = 0; g < 3; ++g ) {
			genotype_levels[g] = double( g ) - mean_genotype ;
		}
	}
	
	using integration::derivative ;
	snptest2::NullModelLogLikelihood null_loglikelihood( phenotype_values ) ;
	integration::Derivative< snptest2::NullModelLogLikelihood > null_loglikelihood_derivative= derivative( null_loglikelihood ) ;
	Vector null_parameters = integration::find_root_by_newton_raphson(
		null_loglikelihood_derivative,
		Vector::Zero( 1 )
	) ;

	m_ui_context.logger() << "AssociationTester:        null: MLE is " << null_parameters << ".\n" ;
	m_ui_context.logger() << "AssociationTester:        null: loglikelihood is " << null_loglikelihood.get_value_of_function() << ".\n" ;

	m_ui_context.logger() << "AssociationTester: alternative: Finding MLE...\n" ;

	snptest2::AlternativeModelLogLikelihood alternative_loglikelihood(
		phenotype_values,
		genotypes,
		genotype_levels
	) ;
	integration::Derivative< snptest2::AlternativeModelLogLikelihood > alternative_loglikelihood_derivative = derivative( alternative_loglikelihood ) ;
	Vector alternative_parameters = integration::find_root_by_newton_raphson(
		alternative_loglikelihood_derivative,
		Vector::Zero( 2 )
	) ;

	m_ui_context.logger() << "AssociationTester: alternative: MLE is " << alternative_parameters << ".\n" ;
	m_ui_context.logger() << "AssociationTester: alternative: loglikelihood is " << alternative_loglikelihood.get_value_of_function() << ".\n" ;
	m_ui_context.logger() << "AssociationTester: -2LR = " << -2.0 * ( null_loglikelihood.get_value_of_function() - alternative_loglikelihood.get_value_of_function() ) << ".\n" ;

	FrequentistTestResults result ;
	result.test_statistic = result.p_value = result.beta = result.standard_error = std::numeric_limits< double >::quiet_NaN() ;
	result.test_statistic = -2.0 * ( null_loglikelihood.get_value_of_function() - alternative_loglikelihood.get_value_of_function() ) ;
	if( result.test_statistic > 0.0 ) {
		result.p_value = boost::math::cdf(
			boost::math::complement(
				m_chi_squared,
				-2.0 * ( null_loglikelihood.get_value_of_function() - alternative_loglikelihood.get_value_of_function() )
			)
		) ;
	}
	result.beta = alternative_parameters( 1 ) ;
	result.variance_covariance = (-alternative_loglikelihood.get_value_of_second_derivative()).inverse() ;
	result.standard_error = std::sqrt( result.variance_covariance( 1, 1 ) ) ;
	return result ;
}

void AssociationTester::end_processing_snps() {
	m_sink.reset() ;
	m_ui_context.logger() << "AssociationTester: finished processing SNPs.\n" ;
}

