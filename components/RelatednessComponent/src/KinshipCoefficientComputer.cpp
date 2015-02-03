
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <boost/ptr_container/ptr_deque.hpp>
#include <boost/function.hpp>
#include <boost/timer/timer.hpp>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include "unistd.h"
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
#include "worker/FunctionTask.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "components/RelatednessComponent/KinshipCoefficientComputer.hpp"
#include "components/RelatednessComponent/PCAComputer.hpp"
#include "components/RelatednessComponent/mean_centre_genotypes.hpp"

#include "boost/threadpool.hpp"

//#define DEBUG_KINSHIP_COEFFICIENT_COMPUTER 2
#define USING_BOOST_THREADPOOL 1

namespace impl {
	#if HAVE_EIGEN
		template< typename Vector, typename Matrix >
		void accumulate_xxt_using_eigen(
			Vector* data,
			Matrix* result,
			int const begin_sample_i, int const end_sample_i,
			int const begin_sample_j, int const end_sample_j,
			double const scale
		) {
			
			if( begin_sample_i == begin_sample_j && end_sample_i == end_sample_j ) {
				result
					->block( begin_sample_i, begin_sample_j, end_sample_i - begin_sample_i, end_sample_j - begin_sample_j )
					.template selfadjointView< Eigen::Lower >()
					.rankUpdate( data->segment( begin_sample_i, end_sample_i - begin_sample_i ), scale ) ;
			} else {
				result
					->block(  begin_sample_i, begin_sample_j, end_sample_i - begin_sample_i, end_sample_j - begin_sample_j ).noalias()
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
			KinshipCoefficientComputer::Computation::Vector* data,
			KinshipCoefficientComputer::Computation::Matrix* result,
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
					end_sample_i - begin_sample_i,
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

	std::ostream& operator<<( std::ostream& out, SampleBounds const& b ) {
		out << "[" << b.begin_sample_i << "-" << b.end_sample_i << ", " << b.begin_sample_j << "-" << b.end_sample_j << "]" ;
		return out ;
	}

	struct ComputeXXtTask: public boost::noncopyable
	{
		typedef KinshipCoefficientComputer::Computation::Vector Vector ;
		typedef KinshipCoefficientComputer::Computation::Matrix Matrix ;
		typedef KinshipCoefficientComputer::Computation::IntegerVector IntegerVector ;
		typedef KinshipCoefficientComputer::Computation::IntegerMatrix IntegerMatrix ;
		typedef std::auto_ptr< ComputeXXtTask > UniquePtr ;
		
		virtual ~ComputeXXtTask() {}

		virtual void add_data(
			KinshipCoefficientComputer::Computation::Vector&,
			KinshipCoefficientComputer::Computation::IntegerVector&
		) = 0 ;

		virtual void operator()() = 0 ;
		
		static UniquePtr create(
			Matrix* result,
			IntegerMatrix* counts,
			SampleBounds const bounds,
			std::string const& method
		) ;
	} ;

	struct ComputeXXtTaskImpl: public ComputeXXtTask {
		typedef boost::function< void (
			Vector* data,
			Matrix* result,
			int const begin_sample_i, int const end_sample_i,
			int const begin_sample_j, int const end_sample_j,
			double const scale
		) >
		AccumulateXXt ;

		typedef boost::function< void (
			IntegerVector* data,
			IntegerMatrix* result,
			int const begin_sample_i, int const end_sample_i,
			int const begin_sample_j, int const end_sample_j,
			double const scale
		) >
		AccumulateXXtInteger ;
		
		ComputeXXtTaskImpl(
			Matrix* result,
			IntegerMatrix* counts,
			SampleBounds const bounds,
			AccumulateXXt accumulate_xxt,
			AccumulateXXtInteger accumulate_xxt_integer
		):
			m_result( result ),
			m_counts( counts ),
			m_bounds( bounds ),
			m_accumulate_xxt( accumulate_xxt ),
			m_accumulate_xxt_integer( accumulate_xxt_integer )
		{
			assert( result ) ;
			assert( m_accumulate_xxt ) ;
		}

		void add_data(
			Vector& data,
			IntegerVector& nonmissingness
		) {
			m_data.push_back( &data ) ;
			m_nonmissingness.push_back( &nonmissingness ) ;
		}

		std::size_t number_of_snps() const { return m_data.size() ; }

		void operator()() {
			for( std::size_t snp_i = 0; snp_i < m_data.size(); ++snp_i ) {
				m_accumulate_xxt(
					m_data[ snp_i ],
					m_result,
					m_bounds.begin_sample_i, m_bounds.end_sample_i,
					m_bounds.begin_sample_j, m_bounds.end_sample_j,
					1.0
				) ;

				m_accumulate_xxt_integer(
					m_nonmissingness[ snp_i ],
					m_counts,
					m_bounds.begin_sample_i, m_bounds.end_sample_i,
					m_bounds.begin_sample_j, m_bounds.end_sample_j,
					1.0
				) ;
			}
		}
		
	protected:
		Matrix* m_result ;
		IntegerMatrix* m_counts ;
		SampleBounds const m_bounds ;
		std::vector< Vector* > m_data ;
		std::vector< IntegerVector* > m_nonmissingness ;
		AccumulateXXt m_accumulate_xxt ;
		AccumulateXXtInteger m_accumulate_xxt_integer ;
	} ;
	
	struct ComputeXXtUsingEigenTask: public ComputeXXtTaskImpl {
		ComputeXXtUsingEigenTask(
			Matrix* result,
			IntegerMatrix* counts,
			SampleBounds const bounds
		):
			ComputeXXtTaskImpl(
				result, counts, bounds,
				&impl::accumulate_xxt_using_eigen< Vector, Matrix >, &impl::accumulate_xxt_using_eigen< IntegerVector, IntegerMatrix >
			)
		{}
	} ;

#if HAVE_CBLAS
	struct ComputeXXtUsingCBlasTask: public ComputeXXtTaskImpl {
		ComputeXXtUsingCBlasTask(
			Matrix* result,
			IntegerMatrix* counts,
			SampleBounds const bounds
		):
			ComputeXXtTaskImpl(
				result, counts, bounds,
				&impl::accumulate_xxt_using_cblas, &impl::accumulate_xxt_using_eigen< IntegerVector, IntegerMatrix >
			)
		{}
	} ;
#endif
	
	ComputeXXtTask::UniquePtr ComputeXXtTask::create(
		Matrix* result,
		IntegerMatrix* counts,
		SampleBounds const bounds,
		std::string const& method
	) {
		if( method == "cblas" ) {
			return ComputeXXtTask::UniquePtr(
				new ComputeXXtUsingCBlasTask(
					result, counts, bounds
				)
			) ;
		} else if( method == "eigen" ) {
			return ComputeXXtTask::UniquePtr(
				new ComputeXXtUsingEigenTask(
					result, counts, bounds
				)
			) ;
		} else {
			assert(0) ;
		}
	}

	struct Dispatcher {
		typedef std::auto_ptr< Dispatcher > UniquePtr ;
		typedef KinshipCoefficientComputer::Computation::Vector Vector ;
		typedef KinshipCoefficientComputer::Computation::Matrix Matrix ;
		typedef KinshipCoefficientComputer::Computation::IntegerVector IntegerVector ;
		typedef KinshipCoefficientComputer::Computation::IntegerMatrix IntegerMatrix ;
		typedef boost::function<
			ComputeXXtTask::UniquePtr ( Matrix*, IntegerMatrix*, SampleBounds const )
		> TaskFactory ;

		virtual ~Dispatcher() {}
		virtual void setup( std::size_t number_of_samples ) = 0 ;
		virtual void get_storage( Vector* data, IntegerVector* nonmissingness ) = 0 ;
		virtual void add_data( Vector* data, IntegerVector* nonmissingness ) = 0 ;
		virtual void wait_until_complete() = 0 ;
	} ;



	//
	// Return a tiling of a dxd matrix into (roughly) B blocks.
	// Blocks should be roughly evenly sized.
	// This implementation makes off-diagonal blocks half the size of diagonal blocks.
	std::vector< SampleBounds > get_matrix_lower_diagonal_tiling( std::size_t d, std::size_t B ) {
		// Because off-diagonal blocks take about twice as long, 
		// we make them half the size vertically.
		// If there are (N ( N+1)/2 ) square blocks this makes
		// (N (N+1) ) - N = N^2 blocks in total.
		// So we want to choose N so that N^2 is about twice the number of threads.
		std::size_t const N = B ;
		std::vector< SampleBounds > result ;
		double K = double( d ) / N ;
//#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
		std::cerr << "get_matrix_lower_diagonal_tiling: " << d << "x" << d << " matrix, "
			<< N << " blocks across top row, "
			<< N*(N+1)/2 << " blocks on lower diagonal.\n" ;
//#endif
		// Set up some initial tasks.
		// Do j (the column) as the outer loop
		// so that tasks go in column-major direction.
		for( std::size_t j = 0; j < N; ++j ) {
			for( std::size_t i = 0; i < N; ++i ) {
				if( i == j ) {
					SampleBounds bounds ;
					bounds.begin_sample_i = std::ceil(i*K) ;
					bounds.end_sample_i = std::min( std::size_t( std::ceil( (i+1)*K ) ), d ) ;
					bounds.begin_sample_j = bounds.begin_sample_i ;
					bounds.end_sample_j = bounds.end_sample_i ;
					result.push_back( bounds ) ;
				}
				else {
					// Do two vertical stripes
					SampleBounds bounds ;
					bounds.begin_sample_i = std::ceil(i*K) ;
					bounds.end_sample_i = std::min( std::size_t( std::ceil( (i+1)*K ) ), d ) ;
					bounds.begin_sample_j = std::ceil(j*K) ;
					bounds.end_sample_j = std::ceil((j+0.5)*K) ;
					result.push_back( bounds ) ;

					bounds.begin_sample_j = bounds.end_sample_j ;
					bounds.end_sample_j = std::min( std::size_t( std::ceil( (j+1)*K ) ), d ) ;
					result.push_back( bounds ) ;
				}
			}
		}
		return result ;
	}


	struct MultiThreadedDispatcher: public Dispatcher {
		typedef void (*accumulate_xxt_t)( Vector*, Matrix*, int const begin_sample_i, int const end_sample_i, int const begin_sample_j, int const end_sample_j, double const ) ;
		typedef void (*accumulate_xxt_integer_t)( IntegerVector*, IntegerMatrix*, int const begin_sample_i, int const end_sample_i, int const begin_sample_j, int const end_sample_j, double const ) ;

		MultiThreadedDispatcher(
			Matrix& result,
			IntegerMatrix& nonmissingness,
			worker::Worker* worker,
			TaskFactory task_factory
		):
			m_result( result ),
			m_nonmissingness( nonmissingness ),
			m_worker( worker ),
			m_task_factory( task_factory ),
			m_storage( 100, std::make_pair( Vector(), IntegerVector() ) ),
			m_storage_index( 0 ),
			m_data_index( 0 )
#if USING_BOOST_THREADPOOL
			,m_pool( m_worker->get_number_of_worker_threads() )
#endif
		{
		}

		~MultiThreadedDispatcher() {
		}
		
		void setup( std::size_t number_of_samples ) {
			// If there are N worker threads go for about two jobs per thread.
			std::size_t N = 2 * std::sqrt( m_worker->get_number_of_worker_threads() ) ;

			m_bounds = get_matrix_lower_diagonal_tiling( number_of_samples, N ) ;
			for( std::size_t i = 0; i < m_bounds.size(); ++i ) {
				m_tasks.push_back( 0 ) ;
			}
		}

		void get_storage( Vector* data, IntegerVector* nonmissingness ) {
			data->swap( m_storage[ m_storage_index % m_storage.size() ].first ) ;
			nonmissingness->swap( m_storage[ m_storage_index % m_storage.size() ].second ) ;
		}

		void add_data( Vector* data, IntegerVector* nonmissingness ) {
			// Take the data (without copying it).
			data->swap( m_storage[ m_storage_index % m_storage.size() ].first ) ;
			nonmissingness->swap( m_storage[ m_storage_index % m_storage.size() ].second ) ;
	
			// Add 20 SNPs at a time before submitting.
			bool const schedule = ( m_data_index % 20 ) == 0 ;
			++m_data_index ;

			for( std::size_t task_index = 0; task_index < m_bounds.size(); ++task_index ) {
				std::size_t const bound_index = task_index ;
				
				if( m_tasks.is_null( task_index ) ) {
					ComputeXXtTask::UniquePtr new_task = m_task_factory(
						&m_result,
						&m_nonmissingness,
						m_bounds[ bound_index ]
					) ;
					m_tasks.replace( task_index, new_task ) ;
				}
				m_tasks[ task_index ].add_data(
					m_storage[ m_storage_index % m_storage.size() ].first,
					m_storage[ m_storage_index % m_storage.size() ].second
				) ;

				if( schedule ) {
					m_pool.schedule( boost::bind( &ComputeXXtTask::operator(), &m_tasks[ task_index ] )) ;
				}
			}
			if( schedule ) {
				m_pool.wait() ;
				for( std::size_t i = 0; i < m_tasks.size(); ++i ) {
					m_tasks.replace( i, 0 ) ;
				}
			}

			++m_storage_index ;
		}
	
		void wait_until_complete() {
		}

	private:
		Matrix& m_result ;
		IntegerMatrix& m_nonmissingness ;
		worker::Worker* m_worker ;
		TaskFactory m_task_factory ;
		std::vector< SampleBounds > m_bounds ;
		boost::ptr_deque< boost::nullable< ComputeXXtTask > > m_tasks ;
		std::vector<
			std::pair<
				KinshipCoefficientComputer::Computation::Vector,
				KinshipCoefficientComputer::Computation::IntegerVector
			>
		> m_storage ;
		std::size_t m_storage_index ;
		std::size_t m_data_index ;
		accumulate_xxt_t m_accumulate_xxt ;
		accumulate_xxt_integer_t m_accumulate_xxt_integer ;
#if USING_BOOST_THREADPOOL
		boost::threadpool::pool m_pool ;
#endif
	} ;

	NormaliseGenotypesAndComputeXXt::UniquePtr NormaliseGenotypesAndComputeXXt::create( worker::Worker* worker, std::string const& method ) {
		return NormaliseGenotypesAndComputeXXt::UniquePtr( new NormaliseGenotypesAndComputeXXt( worker, method) ) ;
	}

	NormaliseGenotypesAndComputeXXt::NormaliseGenotypesAndComputeXXt(
		worker::Worker* worker,
		std::string const& method
	):
		m_call_threshhold( 0.9 ),
		m_allele_frequency_threshhold( 0.001 ),
		m_number_of_snps_included( 0 ),
		m_dispatcher(
			new impl::MultiThreadedDispatcher(
				m_result,
				m_nonmissingness,
				worker,
				boost::bind(
					&ComputeXXtTask::create,
					_1, _2, _3, method
				)
			)
		)
	{}

	KinshipCoefficientComputer::Computation::Matrix const& NormaliseGenotypesAndComputeXXt::result() const { return m_result ; }
	KinshipCoefficientComputer::Computation::IntegerMatrix const& NormaliseGenotypesAndComputeXXt::nonmissingness() const { return m_nonmissingness ; }

	void NormaliseGenotypesAndComputeXXt::begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) {
		m_result.setZero( number_of_samples, number_of_samples ) ;
		m_nonmissingness.setZero( number_of_samples, number_of_samples ) ;
		m_dispatcher->setup( number_of_samples ) ;
	}
	
	void NormaliseGenotypesAndComputeXXt::processed_snp(
		genfile::SNPIdentifyingData const& id_data,
		genfile::VariantDataReader::SharedPtr data_reader
	) {
		KinshipCoefficientComputer::Computation::Vector genotypes ;
		KinshipCoefficientComputer::Computation::IntegerVector nonmissingness ;

		m_dispatcher->get_storage( &genotypes, &nonmissingness ) ;

		genotypes.setZero( m_result.rows() ) ;
		nonmissingness.setZero( m_result.rows() ) ;

		data_reader->get(
			":genotypes:",
			genfile::vcf::get_threshholded_calls( genotypes, nonmissingness, m_call_threshhold, 0, 0, 1, 2 )
		) ;

		double allele_frequency = genotypes.sum() / ( 2.0 * nonmissingness.sum() ) ;
		if( std::min( allele_frequency, 1.0 - allele_frequency ) > m_allele_frequency_threshhold ) {
			pca::mean_centre_genotypes( &genotypes, nonmissingness, allele_frequency ) ;
			// Could compute actual SNP sd:
			//double const sd = std::sqrt( ( genotypes.array().square().sum() - ( genotypes.sum() * genotypes.sum() ) ) / ( nonmissingness.sum() - 1 ) ) ;
			// Or sd based on allele frequency:
			double const sd = std::sqrt( ( 2.0 * allele_frequency * ( 1.0 - allele_frequency )) ) ;
			genotypes /= sd ;
			m_dispatcher->add_data( &genotypes, &nonmissingness ) ;
			++m_number_of_snps_included ;
		}
	}
	
	void NormaliseGenotypesAndComputeXXt::end_processing_snps() {
		m_dispatcher->wait_until_complete() ;
		m_result.array() /= m_nonmissingness.array().cast< double >() ;
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

void KinshipCoefficientComputer::begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& metadata ) {
	m_computation->begin_processing_snps( number_of_samples, metadata ) ;
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
#if 1
	send_results(
		m_computation->nonmissingness().cast< double >(),
		m_computation->result().cast< double >(),
		"KinshipCoefficientComputer",
		description
	) ;
	#endif
}

namespace impl {
	NTaskDispatcher::NTaskDispatcher(
		worker::Worker* worker
	):
		m_worker( worker )
	{
	}
	
	void NTaskDispatcher::set_number_of_tasks( std::size_t number_of_tasks ) {
		wait_until_complete() ;
		// TODO: do this better using the constructor
		// I don't have the docs here so don't know the right call.
		m_tasks.clear() ;
		for( std::size_t i = 0; i < number_of_tasks; ++i ) {
			m_tasks.push_back( 0 ) ;
		}
	}

	std::size_t NTaskDispatcher::number_of_tasks() const { return m_tasks.size() ; }
	
	void NTaskDispatcher::submit_task( std::size_t task_i, boost::function< void() > task ) {
		assert( task_i < m_tasks.size() ) ;
		if( !m_tasks.is_null( task_i ) ) {
			m_tasks[ task_i ].wait_until_complete() ;
		}
		worker::Task::UniquePtr new_task( new worker::FunctionTask( task ) ) ;
		m_tasks.replace( task_i, new_task ) ;
		m_worker->tell_to_perform_task( m_tasks[ task_i ] ) ;
	}
	
	void NTaskDispatcher::wait_until_complete() {
		for( std::size_t i = 0; i < m_tasks.size(); ++i ) {
			if( !m_tasks.is_null( i ) ) {
				m_tasks[ i ].wait_until_complete() ;
			}
		}
	}

	NormaliseGenotypesAndComputeXXtFast::UniquePtr NormaliseGenotypesAndComputeXXtFast::create(
		worker::Worker* worker,
		std::size_t const number_of_snps_per_computation
	) {
		return NormaliseGenotypesAndComputeXXtFast::UniquePtr( new NormaliseGenotypesAndComputeXXtFast( worker, number_of_snps_per_computation ) ) ;
	}

	NormaliseGenotypesAndComputeXXtFast::NormaliseGenotypesAndComputeXXtFast(
		worker::Worker* worker,
		std::size_t const number_of_snps_per_computation
	):
		m_worker( worker ),
		m_dispatcher( new NTaskDispatcher( worker ) ),
		m_call_threshhold( 0.9 ),
		m_allele_frequency_threshhold( 0.001 ),
		m_number_of_snps_per_computation( number_of_snps_per_computation )
	{
		// We allocate lookup exactly once, here.
		// The size must be a power of 4.
		m_lookup_table.resize( std::pow( 4, 2 * m_number_of_snps_per_computation ), 0.0 ) ;
		m_nonmissingness_lookup_table.resize( std::pow( 4, 2 * m_number_of_snps_per_computation ), 1 ) ;
	}

	NormaliseGenotypesAndComputeXXtFast::Matrix const& NormaliseGenotypesAndComputeXXtFast::result() const {
		return m_result ;
	}

	NormaliseGenotypesAndComputeXXtFast::IntegerMatrix const& NormaliseGenotypesAndComputeXXtFast::nonmissingness() const {
		return m_nonmissingness ;
	}

	std::size_t NormaliseGenotypesAndComputeXXtFast::number_of_snps_included() const {
		return m_number_of_snps_included ;
	}

	void NormaliseGenotypesAndComputeXXtFast::begin_processing_snps(
		std::size_t number_of_samples,
		genfile::SNPDataSource::Metadata const&
	) {
		m_result.resize( number_of_samples, number_of_samples ) ;
		m_result.setZero() ;
		m_nonmissingness.resize( number_of_samples, number_of_samples ) ;
		m_nonmissingness.setZero() ;
		std::size_t numberOfTasks = 2 * std::sqrt( m_worker->get_number_of_worker_threads() ) ;
		m_matrix_tiling = get_matrix_lower_diagonal_tiling(
			number_of_samples,
			numberOfTasks
		) ;
		m_dispatcher->set_number_of_tasks( m_matrix_tiling.size() ) ;
		m_combined_genotypes.resize( number_of_samples, 0 ) ;
		m_per_snp_genotypes.resize( number_of_samples, 0 ) ;
		m_number_of_snps_included = 0 ;
	}

	namespace impl {
		template< typename T >
		double add( T a, T b ) {
			return a + b ;
		}
	}

#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
	namespace {
		template< typename T >
		void print_lookup_table( std::ostream& out, std::vector< T > const& lookup_table, std::size_t const max = 64 ) {
			out << "print_lookup_table(): lookup table of size " << lookup_table.size() << ":" ;
			for( std::size_t i = 0; i < std::min( max, lookup_table.size() ); ++i ) {
				if( i % 8 == 0 ) {
					out << "\n" ;
				} else {
					out << " " ;
				}
				out << ( boost::format( " - 0x%02x: %.2f" ) % i % lookup_table[i] ) ;
			}
			out << "\n" ;
		}
	}
#endif

	void NormaliseGenotypesAndComputeXXtFast::processed_snp(
		genfile::SNPIdentifyingData const& id_data,
		genfile::VariantDataReader::SharedPtr data_reader
	) {
		// std::vector< std::size_t > m_genotypes
		// std::vector< int >
		m_per_snp_genotypes.assign( m_result.rows(), 0 ) ;

		// I find it simplest here to encode genotypes as
		// 0 (missing), 1 (AA homozygote), 2 (heterozygote), 3 (BB homozygote).
		data_reader->get(
			":genotypes:",
			genfile::vcf::get_threshholded_calls( m_per_snp_genotypes, m_call_threshhold, 0, 1, 2, 3 )
		) ;

		// compute allele frequency
		double allele_frequency = 0.0 ;
		double nonmissing_count = 0.0 ;
		for( std::size_t i = 0; i < m_combined_genotypes.size(); ++i ) {
			bool const missing = ( m_per_snp_genotypes[i] == 0 ) ;
			if( !missing ) {
				allele_frequency += m_per_snp_genotypes[i] - 1.0 ;
				nonmissing_count += 1 ;
			}
		}
		allele_frequency /= ( nonmissing_count * 2.0 ) ;
		if( std::min( allele_frequency, 1.0 - allele_frequency ) > m_allele_frequency_threshhold ) {
			// Ok we will process this SNP.
			// The three mean-centered genotype levels are 0 (for missing level), -2f, 1-2f, 2-2f
			// We must build a table that handles this.
			// We do this by expanding the table four times (once for each genotype class) on each SNP.
			// 
			std::size_t const current_size = m_lookup_table.size() ;
			std::size_t const lookup_snp_index = m_number_of_snps_included % m_number_of_snps_per_computation ;
			if( lookup_snp_index == 0 ) {
				m_dispatcher->wait_until_complete() ;
				std::copy( m_per_snp_genotypes.begin(), m_per_snp_genotypes.end(), m_combined_genotypes.begin() ) ;
			} else {
				// We allow 4 bits per genotype, although each genotype takes up only the lower two bits.
				// This allows us to encode pairs of genotypes in consecutive pairs of bits
				// for simpler lookup of the lookup table.
				for( std::size_t i = 0; i < m_combined_genotypes.size(); ++i ) {
					m_combined_genotypes[i] |= m_per_snp_genotypes[i] << (4 * lookup_snp_index) ;
				}
			}

			double const mean = 2.0 * allele_frequency ;
			double const sd = std::sqrt( ( 2.0 * allele_frequency * ( 1.0 - allele_frequency )) ) ;
			add_snp_to_lookup_table(
				lookup_snp_index,
				mean, sd,
				&m_lookup_table,
				&m_nonmissingness_lookup_table
			) ;

			++m_number_of_snps_included ;
#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
			boost::format fmt( "%.2f" ) ;
			std::cerr << "NormaliseGenotypesAndComputeXXtFast::processed_snp(): number_of_snps_included: " << m_number_of_snps_included << ".\n" ;
			std::cerr << "NormaliseGenotypesAndComputeXXtFast::processed_snp(): allele_frequency: " << fmt % allele_frequency << ".\n" ;
			std::cerr << "NormaliseGenotypesAndComputeXXtFast::processed_snp(): mean: " << fmt % mean << ", sd = " << fmt % sd << ".\n" ;
			print_lookup_table( std::cerr, m_lookup_table, 256 ) ;
			//print_lookup_table( std::cerr, m_nonmissingness_lookup_table, 256 ) ;
#endif
			if( m_number_of_snps_included % m_number_of_snps_per_computation == 0 ) {
				submit_tasks( m_lookup_table, m_nonmissingness_lookup_table, m_combined_genotypes ) ;
			}
		}
	}

	void NormaliseGenotypesAndComputeXXtFast::add_snp_to_lookup_table(
		std::size_t const snp_lookup_index,
		double const mean,
		double const sd,
		std::vector< double >* lookup_table,
		std::vector< int >* nonmissingness_lookup_table
	) {
		// We lay out the lookup table as follows.
		// The two genotypes at the first SNP occupy the first 4 bits,
		// the two genotypes at the second SNP occupy the next 4 bits,
		// etc.
		// Thus after k SNPs the size of the useful part of the table is 16^k,
		// which is computed as (1 << 4k).
		// Note that because genotype value 0 corresponds to missing data,
		// The block as computed after SNP k-1 is equal to the block
		// for missing values for SNP k.
		if( snp_lookup_index == 0 ) {
			(*nonmissingness_lookup_table)[0x0] = 0 ;
			(*nonmissingness_lookup_table)[0x1] = 0 ;
			(*nonmissingness_lookup_table)[0x2] = 0 ;
			(*nonmissingness_lookup_table)[0x3] = 0 ;
			for( std::size_t g = 1; g < 4; ++g ) {
				(*nonmissingness_lookup_table)[ (g << 2) + 0 ] = 0 ;
				(*lookup_table)[ (g << 2) + 1 ] = ( ( 0.0 - mean ) / sd ) * ( ( (g-1) - mean ) / sd ) ;
				(*lookup_table)[ (g << 2) + 2 ] = ( ( 1.0 - mean ) / sd ) * ( ( (g-1) - mean ) / sd ) ;
				(*lookup_table)[ (g << 2) + 3 ] = ( ( 2.0 - mean ) / sd ) * ( ( (g-1) - mean ) / sd ) ;
			}
		} else {
			// We have already computed the initial block of size 16^k.
			// Now we compute 15 additional blocks for 16 in total.
			// In the computation we break these up as 4 blocks of size 
			// This means this SNP has no effect and we need only copy.
			std::size_t const inner_block_size = 1 << ( 4 * snp_lookup_index ) ;
			std::size_t const outer_block_size = 4 * inner_block_size ;
			for( std::size_t outer_g = 0; outer_g < 4; ++outer_g ) {
				for( std::size_t inner_g = (outer_g == 0) ? 1 : 0; inner_g < 4; ++inner_g ) {
					double const value = ( inner_g == 0 || outer_g == 0 )
						? 0.0
						: (((inner_g-1.0) - mean)/sd) * (((outer_g-1.0) - mean)/sd) ;
					int const nonmissingness = ( inner_g == 0 || outer_g == 0 )
						? 0
						: 1 ;
					std::transform(
						lookup_table->begin(),
						lookup_table->begin() + inner_block_size,
						lookup_table->begin() + outer_block_size * outer_g + inner_block_size * inner_g,
						boost::bind( &impl::add< double >, _1, value )
					) ;
					std::transform(
						nonmissingness_lookup_table->begin(),
						nonmissingness_lookup_table->begin() + inner_block_size,
						nonmissingness_lookup_table->begin() + outer_block_size * outer_g + inner_block_size * inner_g,
						boost::bind( &impl::add< int >, _1, nonmissingness )
					) ;
				}
			}
		}
	}

	void NormaliseGenotypesAndComputeXXtFast::submit_tasks(
		std::vector< double > const& lookup_table,
		std::vector< int > const& nonmissingness_lookup_table,
		std::vector< std::size_t > const& genotypes
	) {
		for( std::size_t tile_i = 0; tile_i < m_matrix_tiling.size(); ++tile_i ) {
	#if 1
			m_dispatcher->submit_task(
				tile_i,
				boost::bind(
					&NormaliseGenotypesAndComputeXXtFast::compute_block,
					this,
					m_matrix_tiling[ tile_i ],
					genotypes,
					m_lookup_table,
					m_nonmissingness_lookup_table
				)
			) ;
	#else
			compute_block( m_matrix_tiling[tile_i], genotypes, lookup_table, nonmissingness_lookup_table ) ;
	#endif
		}
	}

	void NormaliseGenotypesAndComputeXXtFast::end_processing_snps() {
		if( m_number_of_snps_included % m_number_of_snps_per_computation > 0 ) {
			submit_tasks( m_lookup_table, m_nonmissingness_lookup_table, m_combined_genotypes ) ;
		}
		m_dispatcher->wait_until_complete() ;
		m_result.array() /= m_nonmissingness.array().cast< double >() ;
	}

	void NormaliseGenotypesAndComputeXXtFast::compute_block(
		SampleBounds const& sample_bounds,
		std::vector< std::size_t > const& genotypes,
		std::vector< double > const& lookup_table,
		std::vector< int > const& nonmissingness_lookup_table
	) {
#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
		std::cerr << "NormaliseGenotypesAndComputeXXtFast::compute_block(): computing block for " << sample_bounds << ".\n" ;
		std::cerr << "NormaliseGenotypesAndComputeXXtFast::compute_block(): e.g. for samples "
			<< sample_bounds.begin_sample_i << " and " << sample_bounds.begin_sample_j << ":\n"
			<< "lookup_table[" << genotypes[sample_bounds.begin_sample_i] << "] * lookup_table[" << genotypes[sample_bounds.begin_sample_j] << "] = "
			<< lookup_table[ genotypes[sample_bounds.begin_sample_i] ] << " * " << lookup_table[ genotypes[sample_bounds.begin_sample_j] ] << ".\n" ;
				
#endif		
		// only compute the lower diagonal.
		// We also want to make the function go down the rows of the result
		// in the inner loop, since it's stored column-major, so put row index in the inner loop.
		std::size_t const shift = 2 * m_number_of_snps_per_computation ;
		for( std::size_t j = sample_bounds.begin_sample_j; j < sample_bounds.end_sample_j; ++j ) {
			for( std::size_t i = std::max( j, std::size_t( sample_bounds.begin_sample_i ) ); i < sample_bounds.end_sample_i; ++i ) {
				std::size_t const lookup_index = ( genotypes[i] << 2 ) + genotypes[j] ;
				m_result(i,j) += lookup_table[ lookup_index ] ;
				m_nonmissingness(i,j) += nonmissingness_lookup_table[ lookup_index ] ;
			}
		}
	}
}
