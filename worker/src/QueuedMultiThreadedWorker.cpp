
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <sstream>
#include <iomanip>
#include <numeric>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include "worker/QueuedMultiThreadedWorker.hpp"

namespace worker
{
	QueuedMultiThreadedWorker::QueuedMultiThreadedWorker( std::size_t number_of_threads )
		: m_global_state( e_GlobalProcessing ),
		  m_thread_states( number_of_threads, e_Processing ),
		  m_tasks_completed( number_of_threads, 0 ),
		  m_max_queue_size( 100000 )
	{
		for( std::size_t i = 0; i < number_of_threads; ++i ) {
			m_thread_group.create_thread(
				boost::bind(
					&QueuedMultiThreadedWorker::worker_thread_loop,
					this,
					i
				)
			) ;
		}
	}

	QueuedMultiThreadedWorker::~QueuedMultiThreadedWorker() {
		{
			ScopedLock lock( m_queue_mutex ) ;
			m_global_state = e_GlobalTerminating ;
			m_have_task_condition.notify_all() ;
		}
		m_thread_group.join_all() ;		
	}

	bool QueuedMultiThreadedWorker::ask_to_perform_task( Task& task ) {
		ScopedLock lock( m_queue_mutex ) ;
		if( m_task_queue.size() < m_max_queue_size ) {
			m_task_queue.push( &task ) ;
			m_have_task_condition.notify_one() ;
			return true ;
		}
		else {
			return false ;
		}
	}
	
	void QueuedMultiThreadedWorker::tell_to_perform_task( Task& task ) {
		ScopedLock lock( m_queue_mutex ) ;
		while( m_task_queue.size() == m_max_queue_size ) {
			std::cerr << "Waited." ;
			boost::this_thread::sleep( boost::posix_time::milliseconds( 100 ) ) ;
			m_have_capacity_condition.wait( lock ) ;
		}
		m_task_queue.push( &task ) ;
		m_have_task_condition.notify_one() ;
	}

	std::size_t QueuedMultiThreadedWorker::get_number_of_worker_threads() const {
		return m_thread_group.size() ;
	}

	std::string QueuedMultiThreadedWorker::get_summary_of_work_so_far() const {
		std::ostringstream oStream ;
		oStream << "Tasks carried out by " << get_number_of_worker_threads() << " worker threads:\n" ;
		for( std::size_t i = 0; i < get_number_of_worker_threads(); ++i ) {
			oStream << "thread " << std::setw(3) << i << ": " << m_tasks_completed[i] << "\n" ;
		}
		return oStream.str() ;
	}

	std::size_t QueuedMultiThreadedWorker::get_number_of_tasks_completed() const {
		return std::accumulate( m_tasks_completed.begin(), m_tasks_completed.end(), 0 ) ;
	}

	void QueuedMultiThreadedWorker::worker_thread_loop( std::size_t thread_index ) {
		Task* this_thread_task ;
		while( 1 ) {
			if( !( this_thread_task = worker_thread_wait_for_next_task( thread_index ) )) {
				break ;
			}
			worker_thread_perform_task( thread_index, this_thread_task ) ;
		}
	}
	
	Task* QueuedMultiThreadedWorker::worker_thread_wait_for_next_task( std::size_t thread_index ) {
		ScopedLock lock( m_queue_mutex ) ;
		while( m_global_state != e_GlobalTerminating && m_task_queue.empty() ) {
			m_have_task_condition.wait( lock ) ;
		}

		if( m_global_state == e_GlobalTerminating ) {
			return 0 ;
		}
		else {
			return worker_thread_get_task_from_queue() ;
		}
	}
	
	Task* QueuedMultiThreadedWorker::worker_thread_get_task_from_queue() {
		Task* result = m_task_queue.front() ;
		m_task_queue.pop() ;
		m_have_capacity_condition.notify_one() ;
		return result ;
	}
	
	void QueuedMultiThreadedWorker::worker_thread_perform_task( std::size_t thread_index, Task* task ) {
		char& this_thread_state = m_thread_states[ thread_index ] ;
		this_thread_state = e_Processing ;
		perform_task( *task ) ;
		++m_tasks_completed[ thread_index ] ;
		this_thread_state = e_NotProcessing ;
	}
}
