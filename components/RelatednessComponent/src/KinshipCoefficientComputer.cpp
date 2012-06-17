
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
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
#include "genfile/vcf/get_set_eigen.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include "worker/Task.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "components/RelatednessComponent/KinshipCoefficientComputer.hpp"
#include "components/RelatednessComponent/PCAComputer.hpp"
#include "components/RelatednessComponent/mean_centre_genotypes.hpp"
#include "components/RelatednessComponent/mean_centre_genotypes.hpp"

#define DEBUG_KINSHIP_COEFFICIENT_COMPUTER 1

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
		m_finalised( false )
	{}

	void KinshipCoefficientComputerTask::add_snp(
		genfile::SNPIdentifyingData const& id_data,
		genfile::VariantDataReader& data_reader
	) {
		m_id_data.push_back( id_data ) ;
		m_genotype_calls.push_back( Eigen::VectorXd( m_number_of_samples )) ;
		m_non_missing_calls.push_back( Eigen::VectorXd( m_number_of_samples )) ;
		
		// Store missing genotypes as a zero in the m_genotype_calls entry, and in the m_non_missing_calls entry too.
		genfile::vcf::ThreshholdingGenotypeSetter< Eigen::VectorXd > setter( m_genotype_calls.back(), m_non_missing_calls.back(), 0.9, 0, 0, 1, 2 ) ;
		data_reader.get( "genotypes", setter ) ;
	}

	void KinshipCoefficientComputerTask::operator()() {
		for( std::size_t snp_i = 0; snp_i < m_id_data.size(); ++snp_i ) {
			Eigen::VectorXd& genotype_calls = m_genotype_calls[ snp_i ] ;
			Eigen::VectorXd& non_missingness_matrix = m_non_missing_calls[ snp_i ] ;
			double allele_frequency = genotype_calls.sum() / ( 2.0 * m_non_missing_calls[ snp_i ].sum() ) ;
			if( allele_frequency > 0.01 ) {
#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
				std::cerr << std::resetiosflags( std::ios::floatfield ) << std::setprecision( 5 ) ;
				std::cerr << "SNP: " << m_id_data[ snp_i ] << ": freq = " << allele_frequency << ", uncentred genotypes are: " << genotype_calls.transpose().head( 20 ) << "...\n" ;
#endif
				pca::mean_centre_genotypes( &genotype_calls, non_missingness_matrix, allele_frequency ) ;
				
#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
				std::cerr << "mean-centred genotypes are: " << genotype_calls.transpose().head( 20 ) << "...\n" ;
				std::cerr << "non-missingness is: " << non_missingness_matrix.transpose().head( 20 ) << "...\n" ;
#endif
#if HAVE_CBLAS
				// CBLAS is faster for this usage.  Don't know why.
				cblas_dsyr(
					CblasColMajor,
					CblasLower,
					m_number_of_samples,
					1.0 / ( 2.0 * allele_frequency * ( 1.0 - allele_frequency )),
					genotype_calls.data(),
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
				m_result->selfadjointView< Eigen::Lower >().rankUpdate( genotype_calls, 1.0 / ( 2.0 * allele_frequency * ( 1.0 - allele_frequency ) ) ) ;
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
	
	std::string description = "Number of SNPs: "
		+ genfile::string_utils::to_string( m_number_of_snps_processed )
		+ "\nNumber of samples: "
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

