
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef WORKER_TASK_HPP
#define WORKER_TASK_HPP

#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>

namespace worker {
	// A Task is a unit of work.
	// It can be queried to see if it is complete.
	// Tasks are noncopyable by design; this is an abstract
	// base class.  To get work done, derive from this class and
	// override operator().
	struct Task
	{
	public:
		Task() ;
		virtual ~Task() {} ;

	protected:
		virtual void operator()() = 0 ;

	public:
		bool check_if_complete() const ;
		void wait_until_complete() const ;

	private:
		bool m_is_complete ;

		typedef boost::mutex Mutex ;
		typedef Mutex::scoped_lock ScopedLock ;
		typedef boost::condition Condition ;
		mutable Mutex m_mutex ;
		mutable Condition m_condition ;
	private:
		// perform_task calls operator().
		// It then calls set_complete()
		void perform() ;
		// set_complete locks the mutex and sets the completeness flag.
		void set_complete() ;
	private:
		friend class Worker ;
	} ;
}

#endif

