
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
	#if HAVE_EIGEN
		void accumulate_xxt_using_eigen( Eigen::VectorXd* data, Eigen::MatrixXd* result, double const scale ) {
			result
				->selfadjointView< Eigen::Lower >()
				.rankUpdate( *data, scale ) ;
		}
	#endif
	#if HAVE_CBLAS
		void accumulate_xxt_using_cblas( Eigen::VectorXd* data, Eigen::MatrixXd* result, double const scale ) {
			assert( result ) ;
			assert( data ) ;
			int const N = data->size() ;
			assert( N == result->cols() ) ;
			cblas_dsyr(
				CblasColMajor,
				CblasLower,
				N,
				scale,
				data->data(),
				1,
				result->data(),
				N
			) ;
		}
	#endif

	struct ComputeXXtTask: public KinshipCoefficientComputerTask {
		ComputeXXtTask(
			Eigen::MatrixXd* result,
			Eigen::MatrixXd* missing_count,
			AccumulateXXt accumulate_xxt
		):
			m_result( result ),
			m_non_missing_count( missing_count ),
			m_call_threshhold( 0.9 ),
			m_allele_frequency_threshhold( 0.01 ),
			m_accumulate_xxt( accumulate_xxt ),
			m_number_of_samples( result->cols() ),
			m_finalised( false )
		{
			assert( result ) ;
			assert( m_accumulate_xxt ) ;
		}

		void finalise() { m_finalised = true ;}		
		bool is_finalised() const { return m_finalised ;}

		std::size_t number_of_snps() const { return m_data.size() ; }

	protected:
		Eigen::MatrixXd* m_result ;
		Eigen::MatrixXd* m_non_missing_count ;
		double const m_call_threshhold ;
		double const m_allele_frequency_threshhold ;
		std::vector< Eigen::VectorXd* > m_data ;
		std::vector< Eigen::VectorXd* > m_working_space ;
		AccumulateXXt m_accumulate_xxt ;
	private:
		std::size_t const m_number_of_samples ;
		bool m_finalised ;
	} ;

	struct NormaliseGenotypesAndComputeXXtTask: public ComputeXXtTask {
		NormaliseGenotypesAndComputeXXtTask(
			Eigen::MatrixXd* result,
			Eigen::MatrixXd* missing_count,
#if HAVE_CBLAS
			AccumulateXXt accumulate_xxt = &accumulate_xxt_using_cblas
#else
			AccumulateXXt accumulate_xxt = &accumulate_xxt_using_eigen
#endif
		):
			ComputeXXtTask( result, missing_count, accumulate_xxt )
		{
		}

		void add_snp(
			genfile::VariantDataReader::SharedPtr data_reader,
			Eigen::VectorXd& genotypes,
			Eigen::VectorXd& non_missing_calls
		) {
			genfile::vcf::ThreshholdingGenotypeSetter< Eigen::VectorXd > setter( genotypes, non_missing_calls, m_call_threshhold, 0, 0, 1, 2 ) ;
			data_reader->get( "genotypes", setter ) ;
			m_data.push_back( &genotypes ) ;
			m_working_space.push_back( &non_missing_calls ) ;
		}
		
		void operator()() {
			for( std::size_t snp_i = 0; snp_i < m_data.size(); ++snp_i ) {
				Eigen::VectorXd& genotype_calls = *m_data[ snp_i ] ;
				Eigen::VectorXd& non_missing_genotypes = *m_working_space[ snp_i ] ;
				double allele_frequency = genotype_calls.sum() / ( 2.0 * non_missing_genotypes.sum() ) ;
				if( allele_frequency > m_allele_frequency_threshhold ) {
	#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
					std::cerr << std::resetiosflags( std::ios::floatfield ) << std::setprecision( 5 ) ;
					std::cerr << "SNP: freq = " << allele_frequency << ", uncentred genotypes are: " << genotype_calls.transpose().head( 20 ) << "...\n" ;
	#endif
					pca::mean_centre_genotypes( &genotype_calls, non_missing_genotypes, allele_frequency ) ;
					m_accumulate_xxt(
						&genotype_calls,
						m_result,
						1.0 / ( 2.0 * allele_frequency * ( 1.0 - allele_frequency ))
					) ;
					m_accumulate_xxt(
						&non_missing_genotypes,
						m_non_missing_count,
						1.0
					) ;
				}
			}
		}
	} ;

	struct ComputeConcordanceTask: public ComputeXXtTask {
		ComputeConcordanceTask(
			Eigen::MatrixXd* result,
			Eigen::MatrixXd* missing_count,
#if HAVE_CBLAS
			AccumulateXXt accumulate_xxt = &accumulate_xxt_using_cblas
#else
			AccumulateXXt accumulate_xxt = &accumulate_xxt_using_eigen
#endif
		):
			ComputeXXtTask( result, missing_count, accumulate_xxt )
		{
		}

		void add_snp(
			genfile::VariantDataReader::SharedPtr data_reader,
			Eigen::VectorXd& genotypes,
			Eigen::VectorXd& working_space
		) {
			genfile::vcf::ThreshholdingGenotypeSetter< Eigen::VectorXd > setter( genotypes, m_call_threshhold, -1, 0, 1, 2 ) ;
			data_reader->get( "genotypes", setter ) ;
			m_data.push_back( &genotypes ) ;
			m_working_space.push_back( &working_space ) ;
		}
		
		void operator()() {
			for( std::size_t snp_i = 0; snp_i < m_data.size(); ++snp_i ) {
				Eigen::VectorXd& genotype_calls = *m_data[ snp_i ] ;
				Eigen::VectorXd& genotype_indicator = *m_working_space[ snp_i ] ;

				genotype_indicator.resize( genotype_calls.size() ) ;
				for( double g = 0; g < 3; ++g ) {
					genotype_indicator.array() = ( genotype_calls.array() == g ).cast< double >() ;
					m_accumulate_xxt( &genotype_indicator, m_result, 1.0 ) ;
				}
				// deal with missingness.
				genotype_indicator.array() = ( genotype_calls.array() != -1 ).cast< double >() ;
				m_accumulate_xxt( &genotype_indicator, m_non_missing_count, 1.0 ) ;
			}
		}
	} ;
	
	struct ComputeUnnormalisedXXtTask: public ComputeXXtTask {
		struct MatrixSetter ;
		
		ComputeUnnormalisedXXtTask(
			Eigen::MatrixXd* result,
			Eigen::MatrixXd* missing_count,
#if HAVE_CBLAS
			AccumulateXXt accumulate_xxt = &accumulate_xxt_using_cblas,
#else
			AccumulateXXt accumulate_xxt = &accumulate_xxt_using_eigen,
#endif
			std::string const& data_field = "XY"
		):
			ComputeXXtTask( result, missing_count, accumulate_xxt ),
			m_data_field( data_field )
		{
		}

		void add_snp(
			genfile::VariantDataReader::SharedPtr data_reader,
			Eigen::VectorXd& data,
			Eigen::VectorXd& non_missing_data
		) {
			VectorFiller setter( data, non_missing_data, 2, 0 ) ;
			data_reader->get( m_data_field, setter ) ;
			m_data.push_back( &data ) ;
			m_working_space.push_back( &non_missing_data ) ;
		}

		void operator()() {
			for( std::size_t snp_i = 0; snp_i < m_data.size(); ++snp_i ) {
				Eigen::VectorXd& data = *m_data[ snp_i ] ;
				Eigen::VectorXd& non_missing_data = *m_working_space[ snp_i ] ;
				m_accumulate_xxt(
					&data,
					m_result,
					1.0
				) ;
				m_accumulate_xxt(
					&non_missing_data,
					m_non_missing_count,
					1.0
				) ;
			}
		}

	private:
		std::string const m_data_field ;
		Eigen::MatrixXd m_tmp_intensities ;
		Eigen::MatrixXd m_tmp_non_missing_data ;
	private:
		struct VectorFiller: public genfile::VariantDataReader::PerSampleSetter
		{
			typedef Eigen::VectorXd Vector ;
			VectorFiller(
				Vector& result,
				Vector& nonmissingness,
				int const entries_per_sample,
				double const missing_value = 0
			):
				m_result( result ),
				m_nonmissingness( nonmissingness ),
				m_entries_per_sample( entries_per_sample ),
				m_missing_value( missing_value )
			{}

			void set_number_of_samples( std::size_t n ) {
				m_number_of_samples = n ;
				m_result.resize( n * m_entries_per_sample ) ;
				m_nonmissingness.resize( n * m_entries_per_sample ) ;
			}
			void set_sample( std::size_t n ) {
				assert( n < m_number_of_samples ) ; 
				m_sample = n ;
				m_entry_i = 0 ;
			}
			void set_number_of_entries( std::size_t n ) {
				if( n != m_entries_per_sample ) {
					throw genfile::BadArgumentError( "impl::ComputeUnnormalisedXXtTask::VectorFiller::set_number_of_entries()", "n=" + genfile::string_utils::to_string( n ) ) ;
				}
			}

			void operator()( genfile::MissingValue const value ) {
				int const index = m_sample * m_entries_per_sample + m_entry_i++ ;
				m_result( index ) = m_missing_value ;
				m_nonmissingness( index ) = 0 ;
			}
			void operator()( double const value ) {
				int const index = m_sample * m_entries_per_sample + m_entry_i++ ;
				m_result( index ) = m_missing_value ;
				m_nonmissingness( index ) = 1 ;
			}

		private:
			Vector& m_result ;
			Vector& m_nonmissingness ;
			int const m_entries_per_sample ;
			double const m_missing_value ;
			std::size_t m_number_of_samples ;
			std::size_t m_sample ;
			std::size_t m_number_of_entries ;
			std::size_t m_entry_i ;
		} ;
		
	} ;
}

impl::KinshipCoefficientComputerTask::UniquePtr KinshipCoefficientComputer::compute_kinship(
	std::size_t const number_of_samples,
	Eigen::MatrixXd* result,
	Eigen::MatrixXd* non_missingness,
	impl::KinshipCoefficientComputerTask::AccumulateXXt accumulate_xxt,
	bool const initialise = false
) {
	assert( result ) ;
	assert( non_missingness ) ;
	if( initialise ) {
		result->setZero( number_of_samples, number_of_samples ) ;
		non_missingness->setZero( number_of_samples, number_of_samples ) ;
	}
	return impl::KinshipCoefficientComputerTask::UniquePtr(
		new impl::NormaliseGenotypesAndComputeXXtTask(
			result,
			non_missingness,
			accumulate_xxt
		)
	) ;
}

impl::KinshipCoefficientComputerTask::UniquePtr KinshipCoefficientComputer::compute_concordance(
	std::size_t const number_of_samples,
	Eigen::MatrixXd* result,
	Eigen::MatrixXd* non_missingness,
	impl::KinshipCoefficientComputerTask::AccumulateXXt accumulate_xxt,
	bool const initialise = false
) {
	assert( result ) ;
	assert( non_missingness ) ;
	if( initialise ) {
		result->setZero( number_of_samples, number_of_samples ) ;
		non_missingness->setZero( number_of_samples, number_of_samples ) ;
	}
	return impl::KinshipCoefficientComputerTask::UniquePtr(
		new impl::ComputeConcordanceTask(
			result,
			non_missingness,
			accumulate_xxt
		)
	) ;
}

impl::KinshipCoefficientComputerTask::UniquePtr KinshipCoefficientComputer::compute_intensity_covariance(
	std::size_t const number_of_samples,
	Eigen::MatrixXd* result,
	Eigen::MatrixXd* non_missingness,
	impl::KinshipCoefficientComputerTask::AccumulateXXt accumulate_xxt,
	bool const initialise = false
) {
	assert( result ) ;
	assert( non_missingness ) ;
	if( initialise ) {
		result->setZero( number_of_samples * 2, number_of_samples * 2 ) ;
		non_missingness->setZero( number_of_samples * 2, number_of_samples * 2 ) ;
	}
	return impl::KinshipCoefficientComputerTask::UniquePtr(
		new impl::ComputeUnnormalisedXXtTask(
			result,
			non_missingness,
			accumulate_xxt
		)
	) ;
}

KinshipCoefficientComputer::KinshipCoefficientComputer(
	appcontext::OptionProcessor const& options,
	genfile::CohortIndividualSource const& samples,
	worker::Worker* worker,
	appcontext::UIContext& ui_context,
	ComputationFactory computation_factory
):
	m_options( options ),
	m_ui_context( ui_context ),
	m_samples( samples ),
	m_worker( worker ),
	m_computation_factory( computation_factory ),
#if HAVE_CBLAS
	m_accumulate_xxt( options.check( "-use-eigen" ) ? &impl::accumulate_xxt_using_eigen : &impl::accumulate_xxt_using_cblas )
#else
	m_accumulate_xxt( &impl::accumulate_xxt_using_eigen )
#endif
{
	assert( m_worker ) ;
	assert( m_computation_factory ) ;
}

void KinshipCoefficientComputer::begin_processing_snps( std::size_t number_of_samples ) {
	assert( m_samples.get_number_of_individuals() == number_of_samples ) ;
	m_number_of_samples = number_of_samples ;
	m_number_of_snps_processed = 0 ;
	m_number_of_tasks = std::min( m_worker->get_number_of_worker_threads() + 1, 7ul ) ;
	m_number_of_snps_per_task = 100 ;

	m_genotypes.resize( m_number_of_tasks * m_number_of_snps_per_task ) ;
	m_genotype_non_missingness.resize( m_number_of_tasks * m_number_of_snps_per_task ) ;

	m_result.resize( m_number_of_tasks ) ;
	m_non_missing_count.resize( m_number_of_tasks ) ;
	for( std::size_t i = 0; i < m_number_of_tasks; ++i ) {
		m_tasks.push_back(
			m_computation_factory(
				m_number_of_samples,
				&m_result[ i ],
				&m_non_missing_count[ i ],
				m_accumulate_xxt,
				true
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
			m_computation_factory(
				m_number_of_samples,
				&m_result[ m_current_task ],
				&m_non_missing_count[ m_current_task ],
				m_accumulate_xxt,
				false
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
		m_non_missing_count[0],
		m_result[0],
		"KinshipCoefficientComputer",
		description
	) ;
}

