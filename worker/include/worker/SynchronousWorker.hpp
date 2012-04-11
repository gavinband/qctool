
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef WORKER_SYNCHRONOUS_WORKER_HPP
#define WORKER_SYNCHRONOUS_WORKER_HPP

#include <string>
#include <boost/thread/thread.hpp>
#include "worker/Worker.hpp"

namespace worker
{
	// class SynchronousWorker.
	// A Worker which does all its work immediately in the calling thread.
	class SynchronousWorker: public Worker
	{
	public:
		SynchronousWorker() ;
		bool ask_to_perform_task( Task& task ) ;
		void tell_to_perform_task( Task& task ) ;
		std::size_t get_number_of_worker_threads() const { return 1u ; }
		std::string get_summary_of_work_so_far() const ;
		std::size_t get_number_of_tasks_completed() const ;
	private:
		std::size_t m_number_of_tasks_completed ;
	} ;
}
#endif
