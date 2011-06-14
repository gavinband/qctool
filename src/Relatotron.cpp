#include <limits>
#include <numeric>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/bind.hpp>

#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/CohortIndividualSource.hpp"

#include "fputils/floating_point_utils.hpp"

#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "appcontext/OstreamTee.hpp"

#include "worker/Worker.hpp"
#include "worker/FunctionTask.hpp"

#include "string_utils/parse_utils.hpp"

#include "Relatotron.hpp"
#include "SampleBySampleComputation.hpp"
#include "RelatednessBayesFactorComputation.hpp"
#include "ConcordanceComputation.hpp"

void Relatotron::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Relatedness options" ) ;
	options[ "-relatedness" ]
		.set_description( "Compute relatedness matrices pairwise for all samples."
		 	" The value is the stub of the filename containing relatednesses to produce." )
		.set_takes_single_value() ;
	options[ "-concordance" ]
		.set_description( "Compute concordance for all samples." )
		.set_takes_single_value() ;
	options[ "-pairwise-non-missing-count" ]
		.set_description( "Compute pairwise non-missing-call count matrices for all samples." )
		.set_takes_single_value() ;
	options[ "-kinship" ]
		.set_description( "Compute kinship coefficients, as in Powell, Visscher and Goddard (2010) or Astle & Balding (2009)." )
		.set_takes_single_value() ;
	options.option_implies_option( "-relatedness", "-s" ) ;
	options.option_implies_option( "-concordance", "-s" ) ;
	options.option_implies_option( "-kinship", "-s" ) ;

	options[ "-pairwise-sample-rows" ]
		.set_description( "Choose ranges of samples to compute relatedness for."
			" This option should be a comma-separated list of ranges of the form a-b, meaning use all rows between a and b inclusive."
			" A single value means a 1-element range; a range of the form a- or -a means all samples from, or up to"
			" the given one, inclusive." )
		.set_takes_single_value()
		.set_default_value( "0-" ) ;
	options[ "-pairwise-sample-columns" ]
		.set_description( "Choose ranges of samples to compute relatedness for."
			" This option should be a comma-separated list of ranges of the form a-b, meaning use all rows between a and b inclusive."
			" A single value means a 1-element range; a range of the form a- or -a means all samples up to and including the given one." )
		.set_takes_single_value()
		.set_default_value( "0-" ) ;
	options[ "-relatedness-epsilon" ]
		.set_description( "Set the probability of genotyping error at a SNP."
			" This is used to make the model tolerant to genotyping errors." )
		.set_takes_single_value()
		.set_default_value( 1/1000.0 ) ;
	options[ "-relatedness-null" ]
		.set_description( "Set the probability of IBD0, IBD1, and IBD2 in the null model for the Bayes factor."
		 	" These must lie between 0 and 1 and sum to 1." )
		.set_takes_values_per_use( 3 )
		.set_maximum_multiplicity( 1 )
		.set_default_value( 1.0 )
		.set_default_value( 0.0 )
		.set_default_value( 0.0 ) ;
	options[ "-relatedness-alternative" ]
		.set_description( "Set the probability of IBD0, IBD1, and IBD2 in the alternative model for the Bayes factor."
		" These must lie between 0 and 1 and sum to 1." )
		.set_takes_values_per_use( 3 )
		.set_maximum_multiplicity( 1 )
		.set_default_value( 0.0 )
		.set_default_value( 0.0 )	
		.set_default_value( 1.0 ) ;
}

Relatotron::Relatotron( appcontext::OptionProcessor const& options, genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context ):
	m_options( options ),
	m_samples( samples ),
 	m_ui_context( ui_context )
{
	construct_computations() ;
}

void Relatotron::construct_computations() {
	if( m_options.check_if_option_was_supplied( "-relatedness" )) {
		std::string name = "relatedness" ;
		m_computations.insert(
			name,
			SampleBySampleComputation::create( name, m_options, m_ui_context )
		) ;
		m_computation_files[ "relatedness" ] = m_options.get_value< std::string >( "-relatedness" ) ;
	}
	if( m_options.check_if_option_was_supplied( "-kinship" )) {
		std::string name = "kinship" ;
		m_computations.insert(
			name,
			SampleBySampleComputation::create( name, m_options, m_ui_context )
		) ;
		m_computation_files[ name ] = m_options.get_value< std::string >( "-kinship" ) ;
	}
	if( m_options.check_if_option_was_supplied( "-concordance" )) {
		std::string name = "concordance" ;
		m_computations.insert(
			name,
			SampleBySampleComputation::create( name, m_options, m_ui_context )
		) ;
		m_computation_files[ name ] = m_options.get_value< std::string >( "-concordance" ) ;
	}
	if( m_options.check_if_option_was_supplied( "-pairwise-non-missing-count" )) {
		std::string name = "pairwise-non-missing-count" ;
		m_computations.insert(
			name,
			SampleBySampleComputation::create( name, m_options, m_ui_context )
		) ;
		m_computation_files[ name ] = m_options.get_value< std::string >( "-pairwise-non-missing-count" ) ;
	}
}

void Relatotron::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	m_number_of_samples = number_of_samples ;
	m_snps.clear() ;
	m_snps.reserve( number_of_snps ) ;
	m_genotypes.clear() ;
	m_genotypes.reserve( number_of_snps ) ;
}

void Relatotron::processed_snp( SNPIdentifyingData const& id_data, SingleSNPGenotypeProbabilities const& genotypes ) {
	assert( genotypes.get_number_of_samples() == m_number_of_samples ) ;
	assert( m_snps.size() == m_genotypes.size() ) ;
	m_snps.push_back( id_data ) ;
	m_genotypes.push_back( genotypes ) ;
}

void Relatotron::end_processing_snps() {
	assert( m_snps.size() == m_genotypes.size() ) ;

	std::size_t estimated_memory_usage = ( m_snps.capacity() * sizeof( SNPIdentifyingData ) )
		+ ( m_genotypes.capacity() * sizeof( SingleSNPGenotypeProbabilities ))
		+ ( m_genotypes.capacity() * m_number_of_samples * 3 * sizeof( double )) ;
	
	m_ui_context.logger()
		<< "Relatotron: finished loading "
		<< m_snps.size()
		<< " SNPs, memory usage is "
		<< std::fixed << std::setprecision( 1 ) << ( estimated_memory_usage / 1000000.0 )
		<< "Mb.\n" ;
		
	for( Computations::iterator i = m_computations.begin(); i != m_computations.end(); ++i ) {
		i->second->prepare( m_snps, m_genotypes ) ;
	}
}

void Relatotron::process( worker::Worker* worker ) {
	std::vector< std::size_t > row_samples = parse_row_spec( m_options.get_value< std::string >( "-pairwise-sample-rows" )) ;
	std::vector< std::size_t > column_samples = parse_row_spec( m_options.get_value< std::string >( "-pairwise-sample-columns" )) ;

	for( Computations::iterator computation_i = m_computations.begin(); computation_i != m_computations.end(); ++computation_i ) {
		// Store the results in a plain matrix.  This is space-inefficient because we only compute the
		// upper diagonal, but it's the easiest thing to do.

		Matrix result = ConstantMatrix( m_number_of_samples, m_number_of_samples, -std::numeric_limits< double >::infinity() ) ;

		if( worker ) {
			process_multithreaded(
				*(computation_i->second),
				&result,
				row_samples,
				column_samples,
				*worker
			) ;
		}
		else {
			process_singlethreaded(
				*(computation_i->second),
				&result,
				row_samples,
				column_samples
			) ;
		}

		m_ui_context.logger() << "Top left of " << computation_i->first << " matrix: [\n" ;
		for( std::size_t i = 0; i < std::min( std::size_t( 6 ), row_samples.size() ); ++i ) {
			for( std::size_t j = 0; j < std::min( std::size_t( 6 ), column_samples.size() ); ++j ) {
				m_ui_context.logger() << std::setprecision( 2 ) << std::setw(8) ;
				if( result( row_samples[i], column_samples[j] ) == -std::numeric_limits< double >::infinity() ) {
					m_ui_context.logger() << "" << " " ;
				}
				else {
					m_ui_context.logger() << result( row_samples[i], column_samples[j] ) << " " ;
				}
			}
			m_ui_context.logger()  << "\n" ;
		}
		m_ui_context.logger() << "]\n" ;
		
		std::string const filename = m_computation_files.find( computation_i->first )->second ;
		
		write_sample_by_sample_matrix(
			result,
			filename,
			computation_i->second->get_summary(),
			row_samples,
			column_samples
		) ;
	}
}

std::vector< std::size_t > Relatotron::parse_row_spec( std::string const& spec ) const {
	std::set< std::size_t > result ;
	std::vector< std::string > specs = string_utils::split_and_strip( spec, "," ) ;
	for( std::size_t i = 0; i < specs.size(); ++i ) {
		std::vector< std::string > this_spec = string_utils::split_and_strip( specs[i], "-" ) ;
		std::size_t a = 0 ;
		std::size_t b = m_number_of_samples - 1 ;
		if( this_spec.size() == 1 ) {
			this_spec.push_back( this_spec[0] ) ;
		}
		if( this_spec.size() != 2 ) {
			throw genfile::BadArgumentError( "Relatotron::parse_row_spec()", "spec \"" + specs[i] + "\" is malformed." ) ;
		}
		if( this_spec[0].size() > 0 ) {
			a = string_utils::to_repr< std::size_t > ( this_spec[0] ) ;
		}
		if( this_spec[1].size() > 0 ) {
			b = string_utils::to_repr< std::size_t > ( this_spec[1] ) ;
		}
		if( a > b ) {
			throw genfile::BadArgumentError( "Relatotron::parse_row_spec()", "spec \"" + specs[i] + "\" is not a valid range." ) ;
		}
		if( a >= m_number_of_samples || b >= m_number_of_samples ) {
			throw genfile::BadArgumentError( "Relatotron::parse_row_spec()", "spec \"" + specs[i] + "\" goes past maximum sample id " + string_utils::to_string( m_number_of_samples - 1 ) + "." ) ;
		}
			for( std::size_t i = a; i <= b; ++i ) {
			result.insert( i ) ;
		}
	}
	return std::vector< std::size_t >( result.begin(), result.end() ) ;
}

void Relatotron::process_multithreaded(
	SampleBySampleComputation& computation,
	Matrix* result,
	std::vector< std::size_t > const& row_samples,
	std::vector< std::size_t > const& column_samples,
	worker::Worker& worker
) {
	boost::ptr_vector< worker::Task > m_tasks ;
	appcontext::UIContext::ProgressContext progress_context = m_ui_context.get_progress_context( "Performing sample-by-sample computation" ) ;
	progress_context.notify_progress( 0, row_samples.size() ) ;
	// Simple scheme: one task per row sample.
	for( std::size_t i = 0; i < row_samples.size(); ++i ) {
		m_tasks.push_back(
			new worker::FunctionTask(
			 	boost::bind(
					&Relatotron::perform_pairwise_computations,
					this,
					boost::ref( computation ),
					result,
					std::vector< std::size_t >( 1, row_samples[i] ),
					boost::ref( column_samples ),
					(appcontext::UIContext::ProgressContext *) ( 0 )
				)
			)
		) ;
		while( !worker.ask_to_perform_task( m_tasks.back() ) ) {
			boost::this_thread::sleep( boost::posix_time::milliseconds( 10 ) ) ;
		} ;
		progress_context.notify_progress( worker.get_number_of_tasks_completed(), m_tasks.size() ) ;
	}

	// Wait for them all to finish...
	for( std::size_t i = 0; i < m_tasks.size(); ++i ) {
		m_tasks[i].wait_until_complete() ;
		progress_context.notify_progress( worker.get_number_of_tasks_completed(), m_tasks.size() ) ;
	}
	progress_context.notify_progress( worker.get_number_of_tasks_completed(), m_tasks.size() ) ;
}

void Relatotron::process_singlethreaded(
	SampleBySampleComputation& computation,
	Matrix* result,
	std::vector< std::size_t > const& row_samples,
	std::vector< std::size_t > const& column_samples
) {
	appcontext::UIContext::ProgressContext progress_context = m_ui_context.get_progress_context( "Performing sample-by-sample computation" ) ;

	this->perform_pairwise_computations(
		computation,
		result,
		row_samples,
		column_samples,
		&progress_context
	) ;
}

void Relatotron::perform_pairwise_computations(
	SampleBySampleComputation& computation,
	Matrix* result,
	std::vector< std::size_t > const& sample1_choice,
	std::vector< std::size_t > const& sample2_choice,
	appcontext::UIContext::ProgressContext* progress_context
) const {
	assert( result->size1() == m_number_of_samples ) ;
	assert( result->size2() == m_number_of_samples ) ;
	if( progress_context ) {
		progress_context->notify_progress( 0, sample1_choice.size() ) ;
	}
	for( std::size_t sample1_i = 0; sample1_i < sample1_choice.size(); ++sample1_i ) {
		std::size_t const sample1 = sample1_choice[ sample1_i ] ;
		for( std::size_t sample2_i = 0; sample2_i < sample2_choice.size(); ++sample2_i ) {
			std::size_t const sample2 = sample2_choice[ sample2_i ] ;
			if( sample2 >= sample1 ) {
				double this_result = computation( sample1, sample2, m_genotypes ) ; ;
				(*result)( sample1, sample2 ) = this_result ;
			}
		}
		if( progress_context ) {
			progress_context->notify_progress( sample1_i + 1, sample1_choice.size() ) ;
		}
	}
}

void Relatotron::write_sample_by_sample_matrix(
	Matrix const& bf_matrix,
	std::string const& filename,
	std::string const& description,
	std::vector< std::size_t > const& row_samples,
	std::vector< std::size_t > const& column_samples
) const {
	appcontext::OUTPUT_FILE_PTR file = appcontext::open_file_for_output( filename ) ;
	(*file) << "# Written by Relatotron, " << appcontext::get_current_time_as_string() << "\n" ;
	(*file) << "# Description: " << description << "\n" ;
	(*file) << "# Number of SNPs: " << m_snps.size() << ".\n" ;
	if( m_genotypes.size() > 0 ) {
		assert( bf_matrix.size1() == m_genotypes[0].size() ) ;
		assert( bf_matrix.size2() == m_genotypes[0].size() ) ;
		assert( bf_matrix.size2() == m_samples.get_number_of_individuals() ) ;
		
		(*file) << "id," ;
		for( std::size_t sample_i = 0; sample_i < column_samples.size(); ++sample_i ) {
			if( sample_i > 0 ) {
				(*file) << "," ;
			}
			(*file) << m_samples.get_entry( column_samples[sample_i], "id_1" ).as< std::string >() ;
		}
		(*file) << "\n" ;

		for( std::size_t sample_i = 0; sample_i < row_samples.size(); ++sample_i ) {
			(*file) << m_samples.get_entry( row_samples[ sample_i ], "id_1" ).as< std::string >() << "," ;
			for( std::size_t sample_j = 0; sample_j < column_samples.size(); ++sample_j ) {
				if( sample_j > 0 ) {
					(*file) << "," ;
				}
				(*file) << bf_matrix( row_samples[ sample_i ], column_samples[ sample_j ] ) ;
			}
			(*file) << "\n" ;
		}
	}
}
