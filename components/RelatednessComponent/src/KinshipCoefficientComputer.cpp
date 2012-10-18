
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

// #define DEBUG_KINSHIP_COEFFICIENT_COMPUTER 1

namespace impl {
	struct NormaliseGenotypesAndComputeXXtTask: public KinshipCoefficientComputerTask {
		NormaliseGenotypesAndComputeXXtTask(
			Eigen::MatrixXd* result,
			Eigen::MatrixXd* missing_count
		) ;

		void add_snp(
			genfile::VariantDataReader::SharedPtr data_reader,
			Eigen::VectorXd& genotypes,
			Eigen::VectorXd& non_missing_calls
		) ;
		
		void finalise() { m_finalised = true ;}		
		bool is_finalised() const { return m_finalised ;}

		std::size_t number_of_snps() const { return m_genotype_calls.size() ; }

		void operator()() ;
		
	private:
		std::size_t const m_number_of_samples ;
		Eigen::MatrixXd* m_result ;
		Eigen::MatrixXd* m_non_missing_count ;
		double const m_call_threshhold ;
		double const m_allele_frequency_threshhold ;
		std::vector< Eigen::VectorXd* > m_genotype_calls ;
		std::vector< Eigen::VectorXd* > m_non_missing_calls ;
		
		bool m_finalised ;
	} ;

	NormaliseGenotypesAndComputeXXtTask::NormaliseGenotypesAndComputeXXtTask(
		Eigen::MatrixXd* result,
		Eigen::MatrixXd* missing_count
	):
		m_number_of_samples( result->cols() ),
		m_result( result ),
		m_non_missing_count( missing_count ),
		m_call_threshhold( 0.9 ),
		m_allele_frequency_threshhold( 0.01 ),
		m_finalised( false )
	{
		assert( result ) ;
	}

	void NormaliseGenotypesAndComputeXXtTask::add_snp(
		genfile::VariantDataReader::SharedPtr data_reader,
		Eigen::VectorXd& genotypes,
		Eigen::VectorXd& non_missing_calls
	) {
		genfile::vcf::ThreshholdingGenotypeSetter< Eigen::VectorXd > setter( genotypes, non_missing_calls, m_call_threshhold, 0, 0, 1, 2 ) ;
		data_reader->get( "genotypes", setter ) ;
		m_genotype_calls.push_back( &genotypes ) ;
		m_non_missing_calls.push_back( &non_missing_calls ) ;
	}

	void NormaliseGenotypesAndComputeXXtTask::operator()() {
		for( std::size_t snp_i = 0; snp_i < m_genotype_calls.size(); ++snp_i ) {
			Eigen::VectorXd& genotype_calls = *m_genotype_calls[ snp_i ] ;
			Eigen::VectorXd& non_missing_genotypes = *m_non_missing_calls[ snp_i ] ;
			double allele_frequency = genotype_calls.sum() / ( 2.0 * non_missing_genotypes.sum() ) ;
			if( allele_frequency > m_allele_frequency_threshhold ) {
#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
				std::cerr << std::resetiosflags( std::ios::floatfield ) << std::setprecision( 5 ) ;
				std::cerr << "SNP: freq = " << allele_frequency << ", uncentred genotypes are: " << genotype_calls.transpose().head( 20 ) << "...\n" ;
#endif
				pca::mean_centre_genotypes( &genotype_calls, non_missing_genotypes, allele_frequency ) ;

#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
				std::cerr << "mean-centred genotypes are: " << genotype_calls.transpose().head( 20 ) << "...\n" ;
				std::cerr << "non-missingness is: " << non_missing_genotypes.transpose().head( 20 ) << "...\n" ;
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
					non_missing_genotypes.data(),
					1,
					m_non_missing_count->data(),
					m_number_of_samples
				) ;
#else
				m_result->selfadjointView< Eigen::Lower >().rankUpdate( genotype_calls, 1.0 / ( 2.0 * allele_frequency * ( 1.0 - allele_frequency ) ) ) ;
				m_non_missing_count->selfadjointView< Eigen::Lower >().rankUpdate( non_missing_genotypes, 1.0 ) ;
#endif
			} else {
#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
				std::cerr << "KinshipCoefficientComputerTask::operator(): omitted SNP with frequency " << std::resetiosflags( std::ios::floatfield ) << allele_frequency << ".\n" ;
#endif
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
{}

void KinshipCoefficientComputer::begin_processing_snps( std::size_t number_of_samples ) {
	assert( m_samples.get_number_of_individuals() == number_of_samples ) ;
	m_number_of_samples = number_of_samples ;
	m_number_of_snps_processed = 0 ;
	m_number_of_tasks = 7 ;//m_worker->get_number_of_worker_threads() + 5 ;
	m_number_of_snps_per_task = 100 ;

	m_result.resize( m_number_of_tasks ) ;
	m_non_missing_count.resize( m_number_of_tasks ) ;
	for( std::size_t i = 0; i < m_result.size(); ++i ) {
		m_result[i].resize( m_number_of_samples, m_number_of_samples ) ;
		m_result[i].setZero() ;
		m_non_missing_count[i].resize( m_number_of_samples, m_number_of_samples ) ;
		m_non_missing_count[i].setZero() ;
	}

	m_genotypes.resize( m_number_of_tasks * m_number_of_snps_per_task ) ;
	m_genotype_non_missingness.resize( m_number_of_tasks * m_number_of_snps_per_task ) ;

	for( std::size_t i = 0; i < m_number_of_tasks; ++i ) {
		m_tasks.push_back(
			new impl::NormaliseGenotypesAndComputeXXtTask(
				&m_result[ i ],
				&m_non_missing_count[ i ]
			)
		) ;
	}
	m_current_task = 0 ;
}

void KinshipCoefficientComputer::processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader::SharedPtr data_reader ) {
 	if( m_tasks[ m_current_task ].is_finalised() ) {
		m_tasks[ m_current_task ].wait_until_complete() ;
		m_tasks.replace(
			m_current_task,
			new impl::NormaliseGenotypesAndComputeXXtTask(
				&m_result[ m_current_task ],
				&m_non_missing_count[ m_current_task ]
			)
		) ;
	}

	std::size_t data_index = ( m_current_task * m_number_of_snps_per_task ) + m_tasks[ m_current_task ].number_of_snps() ;
	m_tasks[ m_current_task ].add_snp(
		data_reader,
		m_genotypes[ data_index ],
		m_genotype_non_missingness[ data_index ]
	) ;

	if( m_tasks[ m_current_task ].number_of_snps() == m_number_of_snps_per_task ) {
		m_tasks[ m_current_task ].finalise() ;
		m_worker->tell_to_perform_task( m_tasks[ m_current_task ] ) ;
		m_current_task = ( m_current_task + 1 ) % m_number_of_tasks ;
	}
	
	++m_number_of_snps_processed ;
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
		m_result[0],
		"KinshipCoefficientComputer",
		description
	) ;
}

