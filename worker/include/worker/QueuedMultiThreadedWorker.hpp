#ifndef WORKER_QUEUED_MULTI_THREADED_WORKER_HPP
#define WORKER_QUEUED_MULTI_THREADED_WORKER_HPP

#include <string>
#include <queue>
#include <boost/thread/thread.hpp>
#include <boost/thread/condition.hpp>
#include <boost/thread/mutex.hpp>
#include "worker/Worker.hpp"

namespace worker
{
	// class QueuedMultiThreadedWorker
	// A Worker which performs all its work in one of a number of "worker" threads.
	// It accepts new work when one of the worker threads is idle, and rejects
	// it when all worker threads are busy.
	class QueuedMultiThreadedWorker: public Worker
	{
	public:
		QueuedMultiThreadedWorker( std::size_t number_of_threads ) ;
		~QueuedMultiThreadedWorker() ;

		bool ask_to_perform_task( Task& task ) ;
		void tell_to_perform_task( Task& task ) ;

		std::size_t get_number_of_worker_threads() const ;

		std::string get_summary_of_work_so_far() const ;
		std::size_t get_number_of_tasks_completed() const ;
		
	private:
		enum ThreadState { e_NotProcessing = 0, e_Processing = 1 } ;
		enum GlobalState { e_GlobalProcessing, e_GlobalTerminating = 2 } ;
		char m_global_state ;
		std::vector< char > m_thread_states ;
		boost::thread_group	m_thread_group ;
		std::vector< std::size_t > m_tasks_completed ;

		typedef boost::mutex Mutex ;
		typedef Mutex::scoped_lock ScopedLock ;
		typedef boost::condition ConditionVariable ;

		mutable Mutex m_queue_mutex ;
		ConditionVariable m_have_task_condition ;
		ConditionVariable m_have_capacity_condition ;
		std::queue< Task* > m_task_queue ;
		std::size_t const m_max_queue_size ;
		
		// To be run in worker threads.
		void worker_thread_loop( std::size_t thread_index ) ;
		Task* worker_thread_wait_for_next_task( std::size_t thread_index ) ;
		Task* worker_thread_get_task_from_queue() ;
		void worker_thread_perform_task( std::size_t thread_index, Task* task ) ;

		// forbid copying, assignment.
		QueuedMultiThreadedWorker( QueuedMultiThreadedWorker const& other ) ;
		QueuedMultiThreadedWorker& operator=( QueuedMultiThreadedWorker const& other ) ;
	} ;
}

#endif
