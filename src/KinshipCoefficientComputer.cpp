#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/function.hpp>
#include "../config.hpp"
#if HAVE_CBLAS
	#include "cblas.h"
#endif
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include "worker/Task.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "KinshipCoefficientComputer.hpp"
#include "PCAComputer.hpp"

namespace impl {
	template< typename Matrix >
	void write_matrix_as_csv(
		std::string const& filename,
		Matrix const& matrix,
		std::string const& source,
		std::string const& description,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_row_names = 0,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_column_names = 0
	) {
		appcontext::OUTPUT_FILE_PTR file = appcontext::open_file_for_output( filename ) ;
		(*file) << "# Created by " << source << ", " << appcontext::get_current_time_as_string() << "\n" ;
		(*file) << description << "\n" ;
		if( get_column_names ) {
			if( get_row_names ) {
				(*file) << "id" ;
			}
			for( int j = 0; j < matrix.cols(); ++j ) {
				if( get_row_names || ( j > 0 )) {
					(*file) << "," ;
				}
				(*file) << get_column_names( j ) ;
			}
			(*file) << "\n" ;
		}

		for( int i = 0; i < matrix.rows(); ++i ) {
			if( (!get_row_names.empty()) ) {
				(*file) << get_row_names(i) ;
			}
			for( int j = 0; j < matrix.cols(); ++j ) {
				if( (!get_row_names.empty()) || ( j > 0 ) ) {
					(*file) << "," ;
				}
				double const& value = matrix( i, j ) ;
				if( value == value ) {
					(*file) << value ;
				}
				else {
					(*file) << "NA" ;
				}
			}
			(*file) << "\n" ;
		}
	}

	template< typename Vector >
	void write_snp_and_vector( boost::shared_ptr< statfile::BuiltInTypeStatSink > sink, std::string const& name, genfile::SNPIdentifyingData snp, Vector const& vector, boost::function< genfile::VariantEntry ( std::size_t ) > get_names ) {
		if( sink->number_of_rows_written() == 0 && sink->current_column() == 0 ) {
			(*sink) | "SNPID" | "rsid" | "chromosome" | "position" | "allele_A" | "allele_B" ;
			if( get_names ) {
				for( int i = 0; i < vector.size(); ++i ) {
					(*sink) | get_names( i ).as< std::string >() ;
				}
			}
			else {
				for( int i = 0; i < vector.size(); ++i ) {
					(*sink) | ( "v" + genfile::string_utils::to_string( i )) ;
				}
			}
		}
		(*sink)
			<< snp.get_SNPID()
			<< snp.get_rsid()
			<< std::string( snp.get_position().chromosome() )
			<< snp.get_position().position()
			<< snp.get_first_allele()
			<< snp.get_second_allele() ;
		for( int i = 0; i < vector.size(); ++i ) {
			(*sink) << vector(i) ;
		}
		(*sink) << statfile::end_row() ;
	}
	
	void write_matrix(
		std::string const& filename,
		Eigen::MatrixXd const& matrix,
		std::string const& source,
		std::string const& description,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_row_names,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_column_names
	) {
		statfile::BuiltInTypeStatSink::UniquePtr sink = statfile::BuiltInTypeStatSink::open( filename ) ;
		sink->write_metadata( "Created by " + source + ", " + appcontext::get_current_time_as_string() ) ;
		sink->write_metadata( "# Description: " + description ) ;

		if( get_row_names ) {
			(*sink) | "id" ;
		}
		for( int j = 0; j < matrix.cols(); ++j ) {
			(*sink) | genfile::string_utils::to_string( get_column_names( j ) ) ;
		}

		for( int i = 0; i < matrix.rows(); ++i ) {
			if( get_row_names ) {
				(*sink) << genfile::string_utils::to_string( get_row_names(i) ) ;
			}
			for( int j = 0; j < matrix.cols(); ++j ) {
				double const& value = matrix( i, j ) ;
				if( value == value ) {
					(*sink) << value ;
				}
				else {
					(*sink) << "NA" ;
				}
			}
			(*sink) << statfile::end_row() ;
		}
	}

	void threshhold_genotypes(
		genfile::SingleSNPGenotypeProbabilities const& genotypes,
		Eigen::VectorXd* threshholded_genotypes,
		Eigen::VectorXd* non_missingness_matrix,
		double* allele_sum,
		double threshhold
	) {
		std::size_t const number_of_samples = genotypes.get_number_of_samples() ;
		threshholded_genotypes->setConstant( number_of_samples, std::numeric_limits< double >::quiet_NaN() ) ;
		non_missingness_matrix->setZero( number_of_samples ) ;
		(*allele_sum) = 0.0 ;
		for( std::size_t sample_i = 0; sample_i < number_of_samples; ++sample_i ) {
			for( std::size_t g = 0; g < 3; ++g ) {
				if( genotypes( sample_i, g ) >= threshhold ) {
					(*threshholded_genotypes)( sample_i ) = double( g ) ;
					(*non_missingness_matrix)( sample_i ) = 1.0 ;
					(*allele_sum) += g ;
					break ;
				}
			}
		}
	}

	void mean_centre_genotypes( 
		Eigen::VectorXd* threshholded_genotypes,
		Eigen::VectorXd const& non_missingness_matrix,
		double allele_frequency
	) {
		std::size_t const number_of_samples = threshholded_genotypes->size() ;
		assert( std::size_t( non_missingness_matrix.size() ) == number_of_samples ) ;
		for( std::size_t sample_i = 0; sample_i < number_of_samples; ++sample_i ) {
			if( non_missingness_matrix( sample_i ) ) {
				(*threshholded_genotypes)( sample_i ) -= 2.0 * allele_frequency ;
			}
			else {
				(*threshholded_genotypes)( sample_i ) = 0.0 ; // this sample does not contribute for this SNP.
			}
		}
	}

	KinshipCoefficientComputerTask::KinshipCoefficientComputerTask(
		std::size_t number_of_samples,
		Eigen::MatrixXd* result,
		Eigen::MatrixXd* missing_count
	):
		m_number_of_samples( number_of_samples ),
		m_result( result ),
		m_non_missing_count( missing_count ),
		m_threshhold( 0.9 ),
		m_data( Eigen::VectorXd::Constant( m_number_of_samples, std::numeric_limits< double >::quiet_NaN() ) ),
		m_non_missingness_matrix( Eigen::VectorXd::Zero( m_number_of_samples ) ),
		m_finalised( false )
	{}

	void KinshipCoefficientComputerTask::add_snp(
		genfile::SNPIdentifyingData const& id_data,
		genfile::VariantDataReader& data_reader
	) {
		m_id_data.push_back( id_data ) ;
		m_genotypes.push_back( genfile::SingleSNPGenotypeProbabilities( m_number_of_samples ) ) ;
		data_reader.get( "genotypes", m_genotypes.back() ) ;
	}

	void KinshipCoefficientComputerTask::operator()() {
		Eigen::VectorXd& data = m_data ;
		Eigen::VectorXd& non_missingness_matrix = m_non_missingness_matrix ;
		for( std::size_t snp_i = 0; snp_i < m_id_data.size(); ++snp_i ) {
			double allele_sum ;
			threshhold_genotypes( m_genotypes[ snp_i ], &m_data, &non_missingness_matrix, &allele_sum, m_threshhold ) ;
			double const allele_frequency = allele_sum / ( 2.0 * non_missingness_matrix.sum() ) ;
			mean_centre_genotypes( &data, non_missingness_matrix, allele_frequency ) ;

	#if HAVE_CBLAS
			// CBLAS is faster for this usage.  Don't know why.
			cblas_dsyr(
				CblasColMajor,
				CblasLower,
				m_number_of_samples,
				1.0 / ( 2.0 * allele_frequency * ( 1.0 - allele_frequency )),
				data.data(),
				1,
				m_result->data(),
				m_number_of_samples
			) ;

			cblas_dsyr(
				CblasColMajor,
				CblasLower,
				m_number_of_samples,
				1.0,
				non_missingness_matrix.data(),
				1,
				m_non_missing_count->data(),
				m_number_of_samples
			) ;
	#else
			m_result->selfadjointView< Eigen::Lower >().rankUpdate( data, 1.0 / ( 2.0 * allele_frequency * ( 1.0 - allele_frequency ) ) ) ;
			m_non_missing_count->selfadjointView< Eigen::Lower >().rankUpdate( non_missingness_matrix, 1.0 ) ;
	#endif
		}
	}
}

void KinshipCoefficientManager::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Kinship options" ) ;
	options[ "-kinship" ]
		.set_description( "Perform kinship computation using threshholded genotype calls and cblas or Eigen libraries." )
		.set_takes_single_value() ;
	options[ "-load-kinship" ]
		.set_description( "Load a previously-computed kinship matrix from the specified file." )
		.set_takes_single_value() ;
	options[ "-PCA" ]
		.set_description( "Perform eigenvector/eigenvalue decomposition of kinship matrix after it is computed." )
		.set_takes_single_value() ;
	options[ "-PCA-exclusions" ]
		.set_description( "Output a list of exclusions based on outliers in the first N PCA components." )
		.set_takes_single_value() ;
	options[ "-loadings" ]
		.set_description( "Compute loadings in addition to PCA." )
		.set_takes_single_value()
		.set_default_value( 0 ) ;

	options.option_implies_option( "-kinship", "-s" ) ;
	options.option_implies_option( "-load-kinship", "-s" ) ;
	options.option_implies_option( "-load-kinship", "-PCA" ) ;
	options.option_implies_option( "-PCA", "-load-kinship" ) ;
	options.option_implies_option( "-loadings", "-PCA" ) ;
	options.option_excludes_option( "-load-kinship", "-kinship" ) ;
}

KinshipCoefficientManager::UniquePtr KinshipCoefficientManager::create(
	appcontext::OptionProcessor const& options,
	genfile::CohortIndividualSource const& samples,
	worker::Worker* worker,
	appcontext::UIContext& ui_context
) {
	KinshipCoefficientManager::UniquePtr result ;
	if( options.check( "-kinship" )) {
		result.reset( new KinshipCoefficientComputer( options, samples, worker, ui_context ) ) ;
	}
	else if( options.check( "-load-kinship" ) ) {
		result.reset( new PCAComputer( options, samples, worker, ui_context ) ) ;
		if( options.check( "-loadings" )) {
			// Need to set up an output location.
			// This should be set up elsewhere, but we do it here for now.
 			boost::shared_ptr< statfile::BuiltInTypeStatSink > loadings_file(
				statfile::BuiltInTypeStatSink::open( options.get< std::string >( "-PCA" ) + ".loadings.csv" ).release()
			) ;
			result->send_per_variant_results_to(
				boost::bind(
					&impl::write_snp_and_vector< Eigen::VectorXd >,
					loadings_file,
					_1, _2, _3, _4
				)
			) ;
		}
	} else {
		assert(0) ;
	}

	result->send_results_to( &impl::write_matrix_as_csv< Eigen::MatrixXd > ) ;
	return result ;
}

void KinshipCoefficientManager::send_results_to( ResultsCallback callback ) {
	m_result_signal.connect( callback ) ;
}

void KinshipCoefficientManager::send_results(
	std::string const& name,
	Eigen::MatrixXd const& matrix,
	std::string const& source,
	std::string const& description,
	GetNames get_row_names,
	GetNames get_column_names
) {
	m_result_signal( name, matrix, source, description, get_row_names, get_column_names ) ;
}

void KinshipCoefficientManager::send_per_variant_results_to( PerVariantResultsCallback callback ) {
	m_per_variant_result_signal.connect( callback ) ;
}

void KinshipCoefficientManager::send_per_variant_results( std::string const& name, genfile::SNPIdentifyingData const& snp, Eigen::VectorXd const& data, GetNames names ) {
	m_per_variant_result_signal( name, snp, data, names ) ;
}

KinshipCoefficientComputer::KinshipCoefficientComputer(
	appcontext::OptionProcessor const& options,
	genfile::CohortIndividualSource const& samples,
	worker::Worker* worker,
	appcontext::UIContext& ui_context
):
	m_options( options ),
	m_ui_context( ui_context ),
	m_samples( samples ),
	m_worker( worker )
{
	assert( m_options.check( "-kinship" )) ;
	m_filename = m_options.get< std::string >( "-kinship" ) ;
	// remove .csv extension if present.
	std::string const extension = ".csv" ;
	if(
		m_filename.size() >= extension.size()
		&& ( m_filename.compare(
			m_filename.size() - extension.size(),
			extension.size(),
			extension
		) == 0 )
	) {
		m_filename.resize( m_filename.size() - extension.size() ) ;
	}
}

void KinshipCoefficientComputer::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	assert( m_samples.get_number_of_individuals() == number_of_samples ) ;
	m_number_of_samples = number_of_samples ;
	m_number_of_snps = number_of_snps ;
	m_number_of_snps_processed = 0 ;
	m_number_of_tasks = 7 ;//m_worker->get_number_of_worker_threads() + 5 ;
	m_number_of_snps_per_task = 10 ;
	m_result.resize( m_number_of_tasks ) ;
	m_non_missing_count.resize( m_number_of_tasks ) ;
	for( std::size_t i = 0; i < m_result.size(); ++i ) {
		m_result[i].resize( m_number_of_samples, m_number_of_samples ) ;
		m_result[i].setZero() ;
		m_non_missing_count[i].resize( m_number_of_samples, m_number_of_samples ) ;
		m_non_missing_count[i].setZero() ;
	}
	for( std::size_t i = 0; i < m_number_of_tasks; ++i ) {
		m_tasks.push_back(
			new impl::KinshipCoefficientComputerTask(
				m_number_of_samples,
				&m_result[ i ],
				&m_non_missing_count[ i ]
			)
		) ;
	}
	m_current_task = 0 ;
#if HAVE_CBLAS
std::cerr << "KinshipComputer: starting processing, using cblas.\n" ;
#else
std::cerr << "KinshipComputer: starting processing, using Eigen.\n" ;
#endif
}

void KinshipCoefficientComputer::processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader& data_reader ) {
	if( m_tasks[ m_current_task ].is_finalised() ) {
		m_tasks[ m_current_task ].wait_until_complete() ;
		m_tasks.replace(
			m_current_task,
			new impl::KinshipCoefficientComputerTask(
				m_number_of_samples,
				&m_result[ m_current_task ],
				&m_non_missing_count[ m_current_task ]
			)
		) ;
	}

	m_tasks[ m_current_task ].add_snp(
		id_data,
		data_reader
	) ;

	if(
		(++m_number_of_snps_processed == m_number_of_snps )
		|| m_tasks[ m_current_task ].number_of_snps() == m_number_of_snps_per_task
	) {
		m_tasks[ m_current_task ].finalise() ;
		m_worker->tell_to_perform_task( m_tasks[ m_current_task ] ) ;
		m_current_task = ( m_current_task + 3 ) % m_number_of_tasks ;
	}
}

namespace impl {
	genfile::VariantEntry get_eigendecomposition_header( std::string const& prefix, int value ) {
		if( value == 0 ) {
			return genfile::VariantEntry( "eigenvalue" ) ;
		}
		else {
			return genfile::VariantEntry( prefix + genfile::string_utils::to_string( value ) ) ;
		}
	}
	
	std::string get_concatenated_ids( genfile::CohortIndividualSource const* samples, std::size_t i ) {
		return samples->get_entry( i, "id_1" ).as< std::string >() + ":" + samples->get_entry( i, "id_2" ).as< std::string >() ;
	}

	std::string string_and_number( std::string const& s, std::size_t i ) {
		return s + genfile::string_utils::to_string( i ) ;
	}
}

void KinshipCoefficientComputer::end_processing_snps() {
	for( std::size_t i = 0; i < m_tasks.size(); ++i ) {
		if( m_tasks[ i ].is_finalised() ) {
			m_tasks[i].wait_until_complete() ;
		}
	}
	for( std::size_t i = 1; i < m_result.size(); ++i ) {
		m_result[0].noalias() += m_result[i] ;
		m_non_missing_count[0].noalias() += m_non_missing_count[i] ;
	}
	m_result[0].array() /= m_non_missing_count[0].array() ;
	
	std::string description = "# Number of SNPs: "
		+ genfile::string_utils::to_string( m_number_of_snps )
		+ "\n# Number of samples: "
		+ genfile::string_utils::to_string( m_number_of_samples ) ;

	send_results(
		m_filename + ".csv",
		m_result[0],
		"KinshipCoefficientComputer",
		description,
		boost::bind(
			&impl::get_concatenated_ids,
			&m_samples,
			_1
		),
		boost::bind(
			&impl::get_concatenated_ids,
			&m_samples,
			_1
		)
	) ;
}


void PCAComputer::load_matrix( std::string const& filename, Eigen::MatrixXd* matrix ) const {
	statfile::BuiltInTypeStatSource::UniquePtr source = statfile::BuiltInTypeStatSource::open( filename ) ;
	matrix->resize( m_samples.get_number_of_individuals(), m_samples.get_number_of_individuals() ) ;
	std::vector< std::size_t > sample_column_indices( m_samples.get_number_of_individuals() ) ;
	for( std::size_t sample_i = 0; sample_i < m_samples.get_number_of_individuals(); ++sample_i ) {
		sample_column_indices[sample_i] = source->index_of_column( impl::get_concatenated_ids( &m_samples, sample_i )) ;
		if( sample_i > 0 ) {
			if( sample_column_indices[sample_i] <= sample_column_indices[sample_i-1] ) {
				m_ui_context.logger() << "!! Error in PCAComputer::load_matrix(): column order does not match samples.\n" ;
				throw genfile::MalformedInputError( source->get_source_spec(), 0, sample_column_indices[ sample_i ] ) ;
			}
		}
	}

	for( std::size_t sample_i = 0; sample_i < m_samples.get_number_of_individuals(); ++sample_i ) {
		std::string id ;
		// find row corresponding to next sample.
		for( (*source) >> id; (*source) && id != impl::get_concatenated_ids( &m_samples, sample_i ); (*source) >> statfile::end_row() >> id ) ;
		if( !(*source) ) {
			throw genfile::MalformedInputError( source->get_source_spec(), source->number_of_rows_read(), 0 ) ;
		}
		for( std::size_t sample_j = 0; sample_j < m_samples.get_number_of_individuals(); ++sample_j ) {
			(*source)
				>> statfile::ignore( sample_column_indices[ sample_j ] - source->current_column() )
				>> (*matrix)( sample_i, sample_j )
			;
		}
		(*source) >> statfile::ignore_all() ;
	}
	assert( source->number_of_rows_read() == m_samples.get_number_of_individuals() ) ;
}


void PCAComputer::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	m_number_of_snps = number_of_snps ;
	m_number_of_samples = number_of_samples ;
	m_number_of_snps_processed = 0 ;

	if( m_options.check_if_option_was_supplied( "-PCA" )) {
		m_ui_context.logger() << "PCAComputer: Computing eigenvalue decomposition of kinship matrix...\n" ;
		m_solver.compute( m_kinship_matrix ) ;
		m_ui_context.logger() << "PCAComputer: Done, writing results...\n" ;
		Eigen::MatrixXd eigendecomposition( m_number_of_samples, m_number_of_samples + 1 ) ;
		eigendecomposition.block( 0, 0, m_number_of_samples, 1 ) = m_solver.eigenvalues().reverse() ;
		eigendecomposition.block( 0, 1, m_number_of_samples, m_number_of_samples ) = Eigen::Reverse< Eigen::MatrixXd, Eigen::Horizontal >( m_solver.eigenvectors() ) ;

		std::string description = "# Number of SNPs: "
			+ genfile::string_utils::to_string( m_number_of_snps )
			+ "\n# Number of samples: "
			+ genfile::string_utils::to_string( m_number_of_samples )
			+ "\n# Note: the first column contains the eigenvalues."
			+ "\n# Note: column (i+1) contains the eigenvector corresponding to the ith eigenvalue." ;

		std::string filename = m_options.get< std::string >( "-PCA" ) ;
		send_results(
			filename + ".eigendecomposition.csv",
			eigendecomposition,
			"PCAComputer",
			description,
			boost::function< genfile::VariantEntry ( int ) >(),
			boost::bind(
				&impl::get_eigendecomposition_header,
				"v",
				_1
			)
		) ;
	}
	if( m_number_of_eigenvectors_to_compute > 0 ) {
		m_number_of_eigenvectors_to_compute = std::min( m_number_of_eigenvectors_to_compute, m_number_of_samples ) ;
		m_number_of_eigenvectors_to_compute = std::min( m_number_of_eigenvectors_to_compute, m_number_of_snps ) ;
		m_eigenvectors.resize( m_number_of_eigenvectors_to_compute ) ;
	}
	if( m_options.check( "-PCA-exclusions" )) {
		
	}
}

void PCAComputer::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	if( m_number_of_eigenvectors_to_compute > 0 ) {
		m_genotype_probabilities.resize( m_number_of_samples ) ;
		data_reader.get( "genotypes", m_genotype_probabilities ) ;
		double allele_sum ;
		impl::threshhold_genotypes( m_genotype_probabilities, &m_genotype_calls, &m_non_missingness, &allele_sum, m_threshhold ) ;
		double const allele_frequency = allele_sum / ( 2.0 * m_non_missingness.sum() ) ;
		impl::mean_centre_genotypes( &m_genotype_calls, m_non_missingness, allele_frequency ) ;
		//
		// If X is the L\times n matrix (L SNPs, n samples) of (mean-centred, scaled) genotypes, we have
		// computed the eigenvalue decomposition X^t X = U D U^t in m_solver.
		// We want to compute the eigenvectors of X X^t instead.
		// These are given by the columns of S^t, or the rows of S, where
		//   S = X U D^{-1/2}
		// At this point we have a single row of X and so can compute a row of S directly.
		m_eigenvectors =
			(
				m_genotype_calls *
				Eigen::Reverse< Eigen::MatrixXd, Eigen::Horizontal >( m_solver.eigenvectors() )
					.leftCols( m_number_of_eigenvectors_to_compute )
			).array() /
			m_solver.eigenvalues()
				.reverse()
				.head( m_number_of_eigenvectors_to_compute )
				.array()
				.sqrt()
			;
		send_per_variant_results(
			"PCA eigenvectors",
			snp,
			m_eigenvectors,
			boost::bind(
				&impl::string_and_number,
				"v_",
				_1
			)
			) ;
	}
}

