
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <deque>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <functional>
#include "metro/concurrency/threadpool.hpp"

// #define DEBUG 1
namespace metro {
	namespace concurrency {
		threadpool::UniquePtr threadpool::create( int number_of_threads ) {
			return threadpool::UniquePtr(
				new threadpool( number_of_threads )
			) ;
		}

		threadpool::threadpool( int number_of_threads ):
			m_running_tasks( 0 ),
			m_stop( false )
		{
			for( int i = 0; i < number_of_threads; ++i ) {
				m_threads.emplace_back(
					std::bind( &threadpool::thread_loop, this )
				) ;
			}
		}
		
		threadpool::~threadpool() {
			// set stop-condition
			std::unique_lock< std::mutex > lock( m_queue_mutex ) ;
			m_stop = true;
			m_task_available.notify_all();
			lock.unlock();

			// all threads terminate, then we're done.
			for ( auto& thread : m_threads ) {
				thread.join() ;
			}
		}
		
		void threadpool::wait()
		{
			std::unique_lock< std::mutex > lock( m_queue_mutex );
			m_task_finished.wait(
				lock,
				[this](){ return m_tasks.empty() && (m_running_tasks == 0); }
			);
		}

		void threadpool::thread_loop() {
			while(1) {
				std::unique_lock<std::mutex> lock( m_queue_mutex ) ;
				m_task_available.wait(
					lock,
					[this](){ return m_stop || !m_tasks.empty(); }
				) ;
				if( !m_tasks.empty() ) {
					++m_running_tasks;

					auto task = m_tasks.front();
					m_tasks.pop_front();

					{
						lock.unlock();
						task();
					}

					lock.lock();
					--m_running_tasks;
					if( m_running_tasks == 0 && m_tasks.size() == 0 ) {
						lock.unlock() ;
						m_task_finished.notify_one() ;
					}
				} else if( m_stop ) {
					break;
				}
			}
		}
	}
}
