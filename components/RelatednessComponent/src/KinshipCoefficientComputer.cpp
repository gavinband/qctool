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
#include "components/RelatednessComponent/KinshipCoefficientComputer.hpp"
#include "components/RelatednessComponent/PCAComputer.hpp"

namespace impl {
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
			
			if( allele_sum > non_missingness_matrix.sum() * 0.01 ) {
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
			} else {
				std::cerr << "KinshipCoefficientComputerTask::operator(): omitted SNP " << m_id_data[ snp_i ]
					<< " with frequency " << std::resetiosflags( std::ios::floatfield ) << allele_frequency << ".\n" ;
			}
		}
	}
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

void KinshipCoefficientComputer::begin_processing_snps( std::size_t number_of_samples ) {
	assert( m_samples.get_number_of_individuals() == number_of_samples ) ;
	m_number_of_samples = number_of_samples ;
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
		m_tasks[ m_current_task ].number_of_snps() == m_number_of_snps_per_task
	) {
		m_tasks[ m_current_task ].finalise() ;
		m_worker->tell_to_perform_task( m_tasks[ m_current_task ] ) ;
		m_current_task = ( m_current_task + 3 ) % m_number_of_tasks ;
	}
	
	++m_number_of_snps_processed ;
}

namespace impl {
	
	std::string get_concatenated_ids( genfile::CohortIndividualSource const* samples, std::size_t i ) {
		return samples->get_entry( i, "id_1" ).as< std::string >() + ":" + samples->get_entry( i, "id_2" ).as< std::string >() ;
	}
}

void KinshipCoefficientComputer::end_processing_snps() {
	if( !m_tasks[ m_current_task ].is_finalised() ) {
		m_tasks[ m_current_task ].finalise() ;
		m_worker->tell_to_perform_task( m_tasks[ m_current_task ] ) ;
		m_current_task = ( m_current_task + 3 ) % m_number_of_tasks ;
	}
	
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
		+ genfile::string_utils::to_string( m_number_of_snps_processed )
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

