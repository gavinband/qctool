
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <sstream>
#include <boost/thread/thread.hpp>
#include "worker/SynchronousWorker.hpp"

namespace worker
{
	SynchronousWorker::SynchronousWorker():
		m_number_of_tasks_completed(0)
	{}
	
	bool SynchronousWorker::ask_to_perform_task( Task& task ) {
		tell_to_perform_task( task ) ;
		return true ;
	}

	void SynchronousWorker::tell_to_perform_task( Task& task ) {
		perform_task( task ) ;
		++m_number_of_tasks_completed ;
	}

	std::string SynchronousWorker::get_summary_of_work_so_far() const {
		std::ostringstream oStream ;
		oStream << "SynchronousWorker: " << m_number_of_tasks_completed << " tasks completed so far." ;
		return oStream.str() ;
	}
	
	std::size_t SynchronousWorker::get_number_of_tasks_completed() const {
		return m_number_of_tasks_completed ;
	}
	
}
