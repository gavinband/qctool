
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <boost/ptr_container/ptr_deque.hpp>
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

//#define DEBUG_KINSHIP_COEFFICIENT_COMPUTER 1

#define BLOCKWISE_PARALLELISM 1
// #define SNPWISE_PARALLELISM 1

namespace impl {
	#if HAVE_EIGEN
		void accumulate_xxt_using_eigen(
			Eigen::VectorXd* data,
			Eigen::MatrixXd* result,
			int const begin_sample_i, int const end_sample_i,
			int const begin_sample_j, int const end_sample_j,
			double const scale
		) {
			if( begin_sample_i == begin_sample_j && end_sample_i == end_sample_j ) {
				result
					->block( begin_sample_i, begin_sample_j, end_sample_i - begin_sample_i, end_sample_j - begin_sample_j )
					.selfadjointView< Eigen::Lower >()
					.rankUpdate( data->segment( begin_sample_i, end_sample_i - begin_sample_i ), scale ) ;
			} else {
				result
					->block(  begin_sample_i, begin_sample_j, end_sample_i - begin_sample_i, end_sample_j - begin_sample_j )
					+= scale * (
						data->segment( begin_sample_i, end_sample_i - begin_sample_i )
						*
						data->segment( begin_sample_j, end_sample_j - begin_sample_j ).transpose()
					) ;
			}
		}
	#endif
	#if HAVE_CBLAS
		void accumulate_xxt_using_cblas(
			Eigen::VectorXd* data,
			Eigen::MatrixXd* result,
			int const begin_sample_i, int const end_sample_i,
			int const begin_sample_j, int const end_sample_j,
			double const scale
		) {
			assert( result ) ;
			assert( data ) ;
			int const N = data->size() ;
			assert( N == result->cols() ) ;
			if( begin_sample_i == begin_sample_j && end_sample_i == end_sample_j ) {
				cblas_dsyr(
					CblasColMajor,
					CblasLower,
					N,
					scale,
					data->segment( begin_sample_i, end_sample_i - begin_sample_i ).data(),
					1,
					result->block( begin_sample_i, begin_sample_j, end_sample_i - begin_sample_i, end_sample_j - begin_sample_j ).data(),
					result->outerStride()
				) ;
			} else {
				cblas_dger(
					CblasColMajor,
					end_sample_i - begin_sample_i,
					end_sample_j - begin_sample_j,
					scale,
					data->segment( begin_sample_i, end_sample_i - begin_sample_i ).data(), 1,
					data->segment( begin_sample_j, end_sample_j - begin_sample_j ).data(), 1,
					result->block( begin_sample_i, begin_sample_j, end_sample_i - begin_sample_i, end_sample_j - begin_sample_j ).data(),
					result->outerStride()
				) ;
			}
		}
	#endif

	struct SampleBounds {
		int begin_sample_i ;
		int end_sample_i ;
		int begin_sample_j ;
		int end_sample_j ;
	} ;
	
	std::ostream& operator<<( std::ostream& out, SampleBounds const& b ) {
		out << "[" << b.begin_sample_i << "-" << b.end_sample_i << ", " << b.begin_sample_j << "-" << b.end_sample_j << "]" ;
		return out ;
	}
	

	struct ComputeXXtTask: public worker::Task {
		typedef std::auto_ptr< ComputeXXtTask > UniquePtr ;
		
		typedef boost::function< void (
			Eigen::VectorXd* data,
			Eigen::MatrixXd* result,
			int const begin_sample_i, int const end_sample_i,
			int const begin_sample_j, int const end_sample_j,
			double const scale
		) >
		AccumulateXXt ;
		
		ComputeXXtTask(
			Eigen::MatrixXd* result,
			SampleBounds const bounds,
			AccumulateXXt accumulate_xxt,
			std::size_t storage_index
		):
			m_result( result ),
			m_bounds( bounds ),
			m_accumulate_xxt( accumulate_xxt ),
			m_storage_index( storage_index )
		{
			assert( result ) ;
			assert( m_accumulate_xxt ) ;
		}

		void add_data(
			Eigen::VectorXd& data
		) {
			m_data.push_back( &data ) ;
		}

		std::size_t number_of_snps() const { return m_data.size() ; }
		std::size_t storage_index() const { return m_storage_index ; }

		void operator()() {
			for( std::size_t snp_i = 0; snp_i < m_data.size(); ++snp_i ) {
#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
				std::cerr << "ComputeXXtTask::operator(): bounds are: " << m_bounds << ".\n" ;
#endif 				
				m_accumulate_xxt(
					m_data[ snp_i ],
					m_result,
					m_bounds.begin_sample_i, m_bounds.end_sample_i,
					m_bounds.begin_sample_j, m_bounds.end_sample_j,
					1.0
				) ;
			}
		}
		
	protected:
		Eigen::MatrixXd* m_result ;
		SampleBounds const m_bounds ;
		std::vector< Eigen::VectorXd* > m_data ;
		AccumulateXXt m_accumulate_xxt ;
		std::size_t const m_storage_index ;
	} ;

	struct Dispatcher {
		typedef std::auto_ptr< Dispatcher > UniquePtr ;
		virtual ~Dispatcher() {}
		virtual void setup( std::size_t number_of_samples ) = 0 ;
		virtual void get_storage( Eigen::VectorXd* data, Eigen::VectorXd* nonmissingness ) = 0 ;
		virtual void add_data( Eigen::VectorXd* data, Eigen::VectorXd* nonmissingness ) = 0 ;
		virtual void wait_until_complete() = 0 ;
	} ;

	struct MultiThreadedDispatcher: public Dispatcher {
		typedef void (*accumulate_xxt_t)( Eigen::VectorXd*, Eigen::MatrixXd*, int const begin_sample_i, int const end_sample_i, int const begin_sample_j, int const end_sample_j, double const ) ;
		
		MultiThreadedDispatcher(
			Eigen::MatrixXd& result,
			Eigen::MatrixXd& nonmissingness,
			worker::Worker* worker,
			accumulate_xxt_t accumulate_xxt = &accumulate_xxt_using_eigen
		):
			m_result( result ),
			m_nonmissingness( nonmissingness ),
			m_worker( worker ),
			m_storage( 100, Eigen::VectorXd() ),
			m_storage_index( 0 ),
			m_accumulate_xxt( accumulate_xxt )
		{}
		
		void setup( std::size_t number_of_samples ) {
			std::size_t B = m_worker->get_number_of_worker_threads() ;
			std::size_t N = std::floor( ( 1.0 + std::sqrt( 1 + 8 * B ) ) / 2 ) ;
			if( m_worker->get_number_of_worker_threads() == 1 ) {
				N = 1 ;
			}
			std::size_t K = std::ceil( double( number_of_samples ) / N ) ;
//#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
			std::cerr << "MultiThreadedDispatcher:setup(): " << number_of_samples << "x" << number_of_samples << " matrix, "
				<< N << " blocks across top row, "
				<< N*(N+1)/2 << " blocks on lower diagonal.\n" ;
//#endif
			// Set up some initial tasks.
			for( std::size_t i = 0; i < N; ++i ) {
				for( std::size_t j = 0; j <= i; ++j ) {
					SampleBounds bounds ;
					bounds.begin_sample_i = (i*K) ;
					bounds.end_sample_i = std::min( (i+1)*K, number_of_samples ) ;
					bounds.begin_sample_j = (j*K) ;
					bounds.end_sample_j = std::min( (j+1)*K, number_of_samples ) ;
					m_bounds.push_back( bounds ) ;
//#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
					std::cerr << "MultiThreadedDispatcher::setup(): added block " << m_bounds.back() << ".\n" ;
//#endif
				}
			}
		}

		void get_storage( Eigen::VectorXd* data, Eigen::VectorXd* nonmissingness ) {
			wait_for_storage( m_storage_index+1 ) ;
			data->swap( m_storage[ m_storage_index % m_storage.size() ] ) ;
			nonmissingness->swap( m_storage[ (m_storage_index+1) % m_storage.size() ] ) ;
		}

		void wait_for_storage( std::size_t storage_index ) {
			// Want to know if there are any tasks still using storage that overlaps with the storage we need.
			while( !m_tasks.empty() && ( m_tasks.front().storage_index() + m_storage.size() <= storage_index ) ) {
				m_tasks.front().wait_until_complete() ;
				m_tasks.pop_front() ;
			}
		}

		void wait_until_complete() {
			while( !m_tasks.empty() ) {
				m_tasks.front().wait_until_complete() ;
				m_tasks.pop_front() ;
			}
		}

		void add_data( Eigen::VectorXd* data, Eigen::VectorXd* nonmissingness ) {
			// Take the data (without copying it).
			data->swap( m_storage[ (m_storage_index) % m_storage.size() ] ) ;
			nonmissingness->swap( m_storage[ (m_storage_index+1) % m_storage.size() ] ) ;
			m_storage_index += 2 ;
			
			for( std::size_t i = 0; i < m_bounds.size(); ++i ) {
				m_tasks.push_back(
					new ComputeXXtTask(
						&m_result,
						m_bounds[i],
						m_accumulate_xxt,
						m_storage_index - 2
					)
				) ;
				m_tasks.back().add_data( m_storage[ (m_storage_index-2) % m_storage.size() ] ) ;
#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
				std::cerr << "MultiThreadedDispatcher:add_data(): dispatching block " << i << " of result to thread pool.\n" ;
#endif

				m_worker->tell_to_perform_task( m_tasks.back() ) ;
				m_tasks.push_back(
					new ComputeXXtTask(
						&m_nonmissingness,
						m_bounds[i],
						m_accumulate_xxt,
						m_storage_index - 1
					)
				) ;
				m_tasks.back().add_data( m_storage[ (m_storage_index-1) % m_storage.size() ] ) ;
				#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
								std::cerr << "MultiThreadedDispatcher:add_data(): dispatching block " << i << " of nonmissingness to thread pool.\n" ;
				#endif
				m_worker->tell_to_perform_task( m_tasks.back() ) ;
			}
		}
	
	private:
		Eigen::MatrixXd& m_result ;
		Eigen::MatrixXd& m_nonmissingness ;
		worker::Worker* m_worker ;
		std::vector< SampleBounds > m_bounds ;
		boost::ptr_deque< ComputeXXtTask > m_tasks ;
		std::vector< Eigen::VectorXd > m_storage ;
		std::size_t m_storage_index ;
		accumulate_xxt_t m_accumulate_xxt ;
	} ;

	NormaliseGenotypesAndComputeXXt::UniquePtr NormaliseGenotypesAndComputeXXt::create( worker::Worker* worker ) {
		return NormaliseGenotypesAndComputeXXt::UniquePtr( new NormaliseGenotypesAndComputeXXt( worker ) ) ;
	}

	NormaliseGenotypesAndComputeXXt::NormaliseGenotypesAndComputeXXt(
		worker::Worker* worker
	):
		m_call_threshhold( 0.9 ),
		m_allele_frequency_threshhold( 0.001 ),
		m_number_of_snps_included( 0 ),
		m_dispatcher( new impl::MultiThreadedDispatcher( m_result, m_nonmissingness, worker ) )
	{}

	Eigen::MatrixXd const& NormaliseGenotypesAndComputeXXt::result() const { return m_result ; }
	Eigen::MatrixXd const& NormaliseGenotypesAndComputeXXt::nonmissingness() const { return m_nonmissingness ; }

	void NormaliseGenotypesAndComputeXXt::begin_processing_snps( std::size_t number_of_samples ) {
		m_result.setZero( number_of_samples, number_of_samples ) ;
		m_nonmissingness.setZero( number_of_samples, number_of_samples ) ;
		m_dispatcher->setup( number_of_samples ) ;
	}
		
	void NormaliseGenotypesAndComputeXXt::processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader::SharedPtr data_reader ) {
		Eigen::VectorXd genotypes ;
		Eigen::VectorXd nonmissingness ;

		m_dispatcher->get_storage( &genotypes, &nonmissingness ) ;

		genotypes.setZero( m_result.rows() ) ;
		nonmissingness.setZero( m_result.rows() ) ;

		{
			genfile::vcf::ThreshholdingGenotypeSetter< Eigen::VectorXd > setter( genotypes, nonmissingness, m_call_threshhold, 0, 0, 1, 2 ) ;
			data_reader->get( "genotypes", setter ) ;
		}

		double allele_frequency = genotypes.sum() / ( 2.0 * nonmissingness.sum() ) ;
		if( std::min( allele_frequency, 1.0 - allele_frequency ) > m_allele_frequency_threshhold ) {
			pca::mean_centre_genotypes( &genotypes, nonmissingness, allele_frequency ) ;
			genotypes /= std::sqrt( ( 2.0 * allele_frequency * ( 1.0 - allele_frequency )) ) ;
			m_dispatcher->add_data( &genotypes, &nonmissingness ) ;
			++m_number_of_snps_included ;
		}
	}
	
	void NormaliseGenotypesAndComputeXXt::end_processing_snps() {
		m_dispatcher->wait_until_complete() ;
		m_result.array() /= m_nonmissingness.array() ;
	}
	
	std::size_t NormaliseGenotypesAndComputeXXt::number_of_snps_included() const {
		return m_number_of_snps_included ;
	}
}

KinshipCoefficientComputer::KinshipCoefficientComputer(
	appcontext::OptionProcessor const& options,
	genfile::CohortIndividualSource const& samples,
	appcontext::UIContext& ui_context,
	Computation::UniquePtr computation
):
	m_options( options ),
	m_ui_context( ui_context ),
	m_samples( samples ),
	m_computation( computation )
{
	assert( m_computation.get() ) ;
}

void KinshipCoefficientComputer::begin_processing_snps( std::size_t number_of_samples ) {
	m_computation->begin_processing_snps( number_of_samples ) ;
}

void KinshipCoefficientComputer::processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader::SharedPtr data_reader ) {
	m_computation->processed_snp( id_data, data_reader ) ;
}

void KinshipCoefficientComputer::end_processing_snps() {
	m_computation->end_processing_snps() ;

	std::string description = "Number of SNPs: "
		+ genfile::string_utils::to_string( m_computation->number_of_snps_included() )
		+ "\nNumber of samples: "
		+ genfile::string_utils::to_string( m_samples.get_number_of_individuals() )
	;

	send_results(
		m_computation->nonmissingness(),
		m_computation->result(),
		"KinshipCoefficientComputer",
		description
	) ;
}

