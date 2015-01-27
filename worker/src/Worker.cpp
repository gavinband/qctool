
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/thread/thread_time.hpp>
#include "worker/Worker.hpp"
#include "worker/SynchronousWorker.hpp"
#include "worker/QueuedMultiThreadedWorker.hpp"

namespace worker {
	Task::Task(): m_is_complete( false ) {}
	
	bool Task::check_if_complete() const {
		ScopedLock lock( m_mutex ) ;
		return m_is_complete ;
	}

	void Task::wait_until_complete() const {
		ScopedLock lock( m_mutex ) ;
		while( !m_is_complete ) {
			m_condition.wait( lock ) ;
		}
	}
	
	void Task::perform() {
		this->operator()() ;
		set_complete() ;
	}
	
	void Task::set_complete() {
		ScopedLock lock( m_mutex ) ;
		m_is_complete = true ;
		m_condition.notify_all() ;
	}

	std::auto_ptr< Worker > Worker::create( unsigned long number_of_threads ) {
		std::auto_ptr< Worker > result ;
		if( number_of_threads == 0 ) {
			result.reset( new worker::SynchronousWorker() ) ;
		}
		else {
			result.reset( new worker::QueuedMultiThreadedWorker( number_of_threads )) ;
		}
		return result ;
	}

	void Worker::perform_task( Task& task ) {
		task.perform() ;
	}
}
