
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_REGRESSION_CONCURRENCY_THREADPOOL_HPP
#define METRO_REGRESSION_CONCURRENCY_THREADPOOL_HPP

#include <vector>
#include <deque>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <functional>

namespace metro {
	namespace concurrency {
		// This threadpool is derived from the implementation
		// by WhozCraig (https://stackoverflow.com/questions/23896421/efficiently-waiting-for-all-tasks-in-a-threadpool-to-finish)
		struct threadpool {
		public:
			typedef std::unique_ptr< threadpool > UniquePtr ;
			static UniquePtr create( int number_of_threads = std::thread::hardware_concurrency() ) ;
			class BatchScheduler ;

		public:
			threadpool( int number_of_threads = std::thread::hardware_concurrency() ) ;
			~threadpool() ;
			
			template<class Task> void schedule( Task&& task ) {
				std::unique_lock< std::mutex > lock( m_queue_mutex ) ;
				m_tasks.emplace_back( std::forward<Task>( task ) );
				m_task_available.notify_one() ;
			}
		
			void wait() ;

			std::size_t number_of_threads() const { return m_threads.size() ; }

			std::unique_ptr< BatchScheduler > batch_scheduler() { return BatchScheduler::UniquePtr( new BatchScheduler( *this )) ; }

		private:
			std::vector< std::thread > m_threads ;
			std::deque< std::function<void()> > m_tasks ;
			std::mutex m_queue_mutex ;
			std::condition_variable m_task_available ;
			std::condition_variable m_task_finished ;
			unsigned int m_running_tasks ;
			bool m_stop ;

		public:
			// class BatchScheduler
			// This allows multiple jobs to be scheduled without notifying threads
			// until the scheduler is destructed
			friend class BatchScheduler ;
			class BatchScheduler {
			public:
				typedef std::unique_ptr< BatchScheduler > UniquePtr ;
			public:
				BatchScheduler( threadpool& pool ):
					m_pool( pool ),
					m_lock( pool.m_queue_mutex )
				{
				}

				~BatchScheduler() {
					m_lock.unlock() ;
					m_pool.m_task_available.notify_all() ;
				}
				
				template<class Task>
				void add( Task&& task ) {
					m_pool.m_tasks.emplace_back( std::forward<Task>( task ) );
				}
				
			private:
				threadpool& m_pool ;
				std::unique_lock< std::mutex > m_lock ;

				// disallow copying
				BatchScheduler( BatchScheduler const& ) ;
				BatchScheduler& operator=( BatchScheduler const& ) ;
			} ;
		
		private:
			void thread_loop() ;
		} ;
	}
}

#endif
