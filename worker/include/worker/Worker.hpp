
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef WORKER_HPP
#define WORKER_HPP

#include <boost/thread/condition.hpp>
#include <iostream>

#include "worker/Task.hpp"

namespace worker
{
	// A class which performs work, as represented by Task objects.
	// (It is expected that all the methods of this class are called from the same thread).
	struct Worker
	{
		typedef std::auto_ptr< Worker > UniquePtr ;
		static UniquePtr create( unsigned long number_of_threads ) ;
		
		virtual ~Worker() {} ;
		// Ask the task performer to perform a task.
		// Return true if the task is accepted, otherwise false.
		// The task object must outlive the task performer.
		virtual bool ask_to_perform_task( Task& task ) = 0 ;
		// Tell the task performer to perform a task.
		// If the worker doesn't have capacity, the calling thread will stall until it does
		// when the task can be handed over.
		virtual void tell_to_perform_task( Task& task ) = 0 ;
		// Return the number of worker threads we have available.
		virtual std::size_t get_number_of_worker_threads() const = 0 ;
		// Return the number of tasks we've completed.
		virtual std::size_t get_number_of_tasks_completed() const = 0 ;
		// Print a human-readable summary of the work carried out so far.
		virtual std::string get_summary_of_work_so_far() const = 0 ;
	protected:
		void perform_task( Task& task ) ;
	} ;
}
#endif
