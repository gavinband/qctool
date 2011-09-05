#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/function.hpp>
#include "../config.hpp"
#if HAVE_CBLAS
	//#include "Accelerate/Accelerate.h"
	#include "cblas.h"
#endif
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "worker/Task.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "KinshipCoefficientComputer.hpp"

namespace impl {
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
			data.setZero() ;
			non_missingness_matrix.setZero() ;

			double allele_sum = 0.0 ;
			for( std::size_t sample_i = 0; sample_i < m_number_of_samples; ++sample_i ) {
				for( std::size_t g = 0; g < 3; ++g ) {
					if( m_genotypes[ snp_i ]( sample_i, g ) >= m_threshhold ) {
						data( sample_i ) = double( g ) ;
						non_missingness_matrix( sample_i ) = 1.0 ;
						allele_sum += g ;
						break ;
					}
				}
			}
			double const allele_freq = allele_sum / ( 2.0 * non_missingness_matrix.sum() ) ;
			for( std::size_t sample_i = 0; sample_i < m_number_of_samples; ++sample_i ) {
				if( non_missingness_matrix( sample_i ) ) {
					data( sample_i ) -= 2.0 * allele_freq ;
				}
				else {
					data( sample_i ) = 0.0 ; // this sample does not contribute for this SNP.
				}
			}
		
	#if HAVE_CBLAS
			// CBLAS is faster for this usage.  Don't know why.
			cblas_dsyr(
				CblasColMajor,
				CblasUpper,
				m_number_of_samples,
				1.0 / ( 2.0 * allele_freq * ( 1.0 - allele_freq )),
				data.data(),
				1,
				m_result->data(),
				m_number_of_samples
			) ;

			cblas_dsyr(
				CblasColMajor,
				CblasUpper,
				m_number_of_samples,
				1.0,
				non_missingness_matrix.data(),
				1,
				m_non_missing_count->data(),
				m_number_of_samples
			) ;
	#else
			m_result->selfadjointView< Eigen::Upper >().rankUpdate( data, 1.0 / ( 2.0 * allele_freq * ( 1.0 - allele_freq ) ) ) ;
			m_non_missing_count->selfadjointView< Eigen::Upper >().rankUpdate( non_missingness_matrix, 1.0 ) ;
	#endif
		}
	}
}

void KinshipCoefficientComputer::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Kinship options" ) ;
	options[ "-kinship" ]
		.set_description( "Perform kinship computation using threshholded genotype calls and cblas or Eigen libraries." )
		.set_takes_single_value() ;
	options[ "-PCA" ]
		.set_description( "Perform eigenvector/eigenvalue decomposition of kinship matrix after it is computed." ) ;

	options.option_implies_option( "-PCA", "-kinship" ) ;
	options.option_implies_option( "-kinship", "-s" ) ;
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
		if( get_row_names.empty() ) std::cerr << "get_row_names is empty.\n" ;
		if( get_column_names.empty() ) std::cerr << "get_column_names is empty.\n" ;
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
	
	genfile::VariantEntry combine_string_and_int( std::string const& prefix, int value ) {
		return genfile::VariantEntry( prefix + genfile::string_utils::to_string( value ) ) ;
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
	
	// currently we only have the upper triangle, fix this.
	for( std::size_t i = 0; i < ( m_number_of_samples - 1 ); ++i ) {
		m_result[0].col(i).tail( m_number_of_samples - i - 1 ) = m_result[0].row(i).tail( m_number_of_samples-i-1 ) ;
		m_non_missing_count[0].col(i).tail( m_number_of_samples-i-1 ) = m_non_missing_count[0].row( i ).tail( m_number_of_samples-i-1 ) ;
	}
	
	std::string description = "# Number of SNPs: "
		+ genfile::string_utils::to_string( m_number_of_snps )
		+ "\n# Number of samples: "
		+ genfile::string_utils::to_string( m_number_of_samples ) ;

	impl::write_matrix_as_csv(
		m_filename + ".csv",
		m_result[0],
		"KinshipCoefficientComputer",
		description,
		boost::bind(
			&genfile::CohortIndividualSource::get_entry,
			&m_samples,
			_1,
			"id_1"
		),
		boost::bind(
			&genfile::CohortIndividualSource::get_entry,
			&m_samples,
			_1,
			"id_1"
		)
	) ;

	if( m_options.check_if_option_was_supplied( "-PCA" )) {
		m_ui_context.logger() << "KinshipCoefficientComputer: Computing eigenvalue decomposition of kinship matrix...\n" ;

		Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > solver( m_result[0] ) ;
		
		m_ui_context.logger() << "KinshipCoefficientComputer: Done, writing results...\n" ;
		
		impl::write_matrix_as_csv(
			m_filename + ".eigenvectors.csv",
			solver.eigenvectors(),
			"KinshipCoefficientComputer",
			description,
			boost::function< genfile::VariantEntry ( int ) >(),
			boost::bind(
				&impl::combine_string_and_int,
				"v",
				_1
			)
		) ;
		
		impl::write_matrix_as_csv(
			m_filename + ".eigenvalues.csv",
			solver.eigenvalues(),
			"KinshipCoefficientComputer",
			description,
			boost::bind(
				&impl::combine_string_and_int,
				"lambda",
				_1
			),
			boost::bind(
				&impl::combine_string_and_int,
				"value",
				_1
			)
		) ;
	}

	m_result_signal( m_result[0], m_non_missing_count[0] ) ;
}

void KinshipCoefficientComputer::write_output() {
}

void KinshipCoefficientComputer::send_results_to( ResultsCallback callback ) {
	m_result_signal.connect( callback ) ;
}	
