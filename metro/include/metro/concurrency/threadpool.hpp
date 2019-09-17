
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
		public:
			threadpool( int number_of_threads = std::thread::hardware_concurrency() ) ;
			~threadpool() ;
			
			template<class Task> void schedule( Task&& task ) {
				std::unique_lock< std::mutex > lock( m_queue_mutex ) ;
				m_tasks.emplace_back( std::forward<Task>( task ) );
				m_task_available.notify_one() ;
			}
		
			void wait() ;

		private:
			std::vector< std::thread > m_threads ;
			std::deque< std::function<void()> > m_tasks ;
			std::mutex m_queue_mutex ;
			std::condition_variable m_task_available ;
			std::condition_variable m_task_finished ;
			unsigned int m_running_tasks ;
			bool m_stop ;
		
		private:
			void thread_loop() ;
		} ;
	}
}

#endif
