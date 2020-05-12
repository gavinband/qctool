
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <memory>
#include <exception>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include "Eigen/Core"
#include "metro/regression/LogLikelihood.hpp"
#include "metro/regression/Design.hpp"
#include "metro/regression/ThreadedLogLikelihood.hpp"
#include "metro/count_range.hpp"
#include "metro/intersect_ranges.hpp"
//#include "boost/threadpool.hpp"

//#define DEBUG 1

namespace metro {
	namespace regression {
		ThreadedLogLikelihood::UniquePtr ThreadedLogLikelihood::create( Design& design, CreateLL create_ll, int number_of_threads ) {
			return ThreadedLogLikelihood::UniquePtr(
				new ThreadedLogLikelihood( design, create_ll, number_of_threads )
			) ;
		}

		ThreadedLogLikelihood::UniquePtr ThreadedLogLikelihood::create( Design::UniquePtr design, CreateLL create_ll, int number_of_threads ) {
			return ThreadedLogLikelihood::UniquePtr(
				new ThreadedLogLikelihood( design, create_ll, number_of_threads )
			) ;
		}

		ThreadedLogLikelihood::UniquePtr ThreadedLogLikelihood::create( Design& design, CreateLL create_ll, metro::concurrency::threadpool& pool ) {
			return ThreadedLogLikelihood::UniquePtr(
				new ThreadedLogLikelihood( design, create_ll, pool )
			) ;
		}

		ThreadedLogLikelihood::UniquePtr ThreadedLogLikelihood::create( Design::UniquePtr design, CreateLL create_ll, metro::concurrency::threadpool& pool ) {
			return ThreadedLogLikelihood::UniquePtr(
				new ThreadedLogLikelihood( design, create_ll, pool )
			) ;
		}

		ThreadedLogLikelihood::~ThreadedLogLikelihood() {
			if( m_pool_owned ) {
				delete m_pool ;
			}
			if( m_design_owned ) {
				delete m_design ;
			}
		}

		ThreadedLogLikelihood::ThreadedLogLikelihood(
			Design::UniquePtr design,
			CreateLL create_ll,
			int number_of_threads
		):
			m_design( design.release() ),
			m_design_owned( true ),
			m_pool( new metro::concurrency::threadpool( number_of_threads )),
			m_pool_owned( true )
		{
			create_lls( create_ll ) ;
		}

		ThreadedLogLikelihood::ThreadedLogLikelihood(
			Design& design,
			CreateLL create_ll,
			int number_of_threads
		):
			m_design( &design ),
			m_design_owned( false ),
			
			m_pool( new metro::concurrency::threadpool( number_of_threads )),
			m_pool_owned( true )
		{
			create_lls( create_ll ) ;
		}

		ThreadedLogLikelihood::ThreadedLogLikelihood(
			Design::UniquePtr design,
			CreateLL create_ll,
			metro::concurrency::threadpool& pool
		):
			m_design( design.release() ),
			m_design_owned( true ),
			m_pool( &pool ),
			m_pool_owned( false )
		{
			create_lls( create_ll ) ;
		}

		ThreadedLogLikelihood::ThreadedLogLikelihood(
			Design& design,
			CreateLL create_ll,
			metro::concurrency::threadpool& pool
		):
			m_design( &design ),
			m_design_owned( false ),
			m_pool( &pool ),
			m_pool_owned( false )
		{
			create_lls( create_ll ) ;
		}
		
		void ThreadedLogLikelihood::create_lls( CreateLL create_ll ) {
			//int const N = metro::impl::count_range( m_design->nonmissing_samples() ) ;
			int const N = m_design->outcome().rows() ;
			std::size_t const threads = m_pool->number_of_threads() ;
			
			// aim to submit jobs for some number of the number of threads,
			// though at least 128 to avoid tiny jobs
			//int const size = int(std::max( 128ul, (N+4*threads-1) / (4*threads))) ;
			//int const size = int(std::max( 128ul, (N+2*threads-1) / (2*threads))) ;
			int const size = int(std::max( 128ul, (N+threads-1) / (threads))) ;

#if DEBUG
			std::cerr << "ThreadedLogLikelihood::create_lls(): N = " << N << " samples...\n" ;
#endif

			if( N == 0 ) {
				// ensure we have at least one ll, even though no samples.
				m_lls.push_back(
					create_ll( *m_design, m_design->nonmissing_samples() )
				) ;
			} else {
#if DEBUG
				std::cerr << "ThreadedLogLikelihood::create_lls(): design nonmissing samples: " << m_design->nonmissing_samples() << ".\n" ;
#endif
				for( int i = 0; i < N; i += size ) {
					//std::vector< metro::SampleRange > ranges = metro::impl::intersect_ranges( m_design->nonmissing_samples(), SampleRange( i, i + size )) ;
					std::vector< metro::SampleRange > ranges( 1, SampleRange( i, std::min( i + size, N ) )) ;
#if DEBUG
					std::cerr << "ThreadedLogLikelihood::create_lls(): creating ll " << m_lls.size() << " with " << metro::impl::count_range( ranges ) << " samples:\n"
						<< ranges << ".\n" ;
#endif
					m_lls.push_back( create_ll( *m_design, ranges )) ;
				}
			}
			
			// sanity check: do parameter names match?
			std::vector< std::string > parameter_names = m_lls[0].get_parameter_names() ;
			for( std::size_t i = 1; i < m_lls.size(); ++i ) {
				if( m_lls[i].get_parameter_names() != parameter_names ) {
					throw std::runtime_error( "parameter names do not match" ) ;
				}
			}
		}

		std::string ThreadedLogLikelihood::get_parameter_name( std::size_t i ) const {
			return m_lls[0].get_parameter_name( i ) ;
		}

		std::vector< std::string > ThreadedLogLikelihood::get_parameter_names() const {
			return m_lls[0].get_parameter_names() ;
		}

		// Return a lx2 matrix identifying the l parameters.
		// For each parameter the row contains:
		// column 0: the index of the outcome column.
		// column 1: the index of the design matrix column.
		int ThreadedLogLikelihood::number_of_parameters() const {
			return m_lls[0].number_of_parameters() ;
		}

		int ThreadedLogLikelihood::number_of_outcomes() const {
			return m_lls[0].number_of_outcomes() ;
		}

		ThreadedLogLikelihood::IntegerMatrix ThreadedLogLikelihood::identify_parameters() const {
			return m_lls[0].identify_parameters() ;
		}

		void ThreadedLogLikelihood::evaluate_at( Point const& parameters, int const numberOfDerivatives ) {
			{
								auto batch = m_pool->batch_scheduler() ;
				for( std::size_t i = 0; i < m_lls.size(); ++i ) {
					batch->add(
//					m_pool->schedule(
						[this,i,parameters,numberOfDerivatives]() {
	//						std::cerr << "RUNNING: " << i << ".\n" ;
							m_lls[i].evaluate_at( parameters, numberOfDerivatives ) ;
	//						std::cerr << "COMPLETED: " << i << ".\n" ;
						}
					) ;
				}
			}
			//std::cerr << "waiting...\n" ;
			m_pool->wait() ;
		}

		void ThreadedLogLikelihood::evaluate( int const numberOfDerivatives ) {
			{
				auto batch = m_pool->batch_scheduler() ;
				for( std::size_t i = 0; i < m_lls.size(); ++i ) {
					batch->add(
//					m_pool->schedule(
						[this,i,numberOfDerivatives]() {
							m_lls[i].evaluate( numberOfDerivatives ) ;
						}
					) ;
				}
			}
//			std::cerr << "waiting...\n" ;
			m_pool->wait() ;
		}

		ThreadedLogLikelihood::Point const& ThreadedLogLikelihood::parameters() const {
			return m_lls[0].parameters() ;
		}
		
		double ThreadedLogLikelihood::get_value_of_function() const {
			double result = m_lls[0].get_value_of_function() ;
			for( std::size_t i = 1; i < m_lls.size(); ++i ) {
				result += m_lls[i].get_value_of_function() ;
			}
			return result ;
		}

		ThreadedLogLikelihood::Vector ThreadedLogLikelihood::get_value_of_first_derivative() const {
			Vector result = m_lls[0].get_value_of_first_derivative() ;
			for( std::size_t i = 1; i < m_lls.size(); ++i ) {
				result += m_lls[i].get_value_of_first_derivative() ;
			}
			return result ;
		}

		ThreadedLogLikelihood::Matrix ThreadedLogLikelihood::get_value_of_second_derivative() const {
			Matrix result = m_lls[0].get_value_of_second_derivative() ;
			for( std::size_t i = 1; i < m_lls.size(); ++i ) {
				result += m_lls[i].get_value_of_second_derivative() ;
			}
			return result ;
		}
		
		std::string ThreadedLogLikelihood::get_summary() const {
			return "ThreadedLogLikelihood:" + m_lls[0].get_summary() ;
		}
	}
}

