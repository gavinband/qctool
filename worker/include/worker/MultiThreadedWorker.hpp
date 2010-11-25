#ifndef WORKER_MULTI_THREADED_WORKER_HPP
#define WORKER_MULTI_THREADED_WORKER_HPP

#include <string>
#include <vector>
#include <boost/thread/thread.hpp>
#include "worker/Worker.hpp"

namespace worker
{
	// class MultiThreadedWorker
	// A Worker which performs all its work in one of a number of "worker" threads.
	// It accepts new work when one of the worker threads is idle, and rejects
	// it when all worker threads are busy.
	class MultiThreadedWorker: public Worker
	{
	public:
		MultiThreadedWorker( std::size_t number_of_threads ) ;
		~MultiThreadedWorker() ;
		
		bool ask_to_perform_task( Task& task ) ;
		void tell_to_perform_task( Task& task ) ;

		std::size_t get_number_of_worker_threads() const ;
		std::size_t get_number_of_tasks_performed_by_thread( std::size_t ) const ;

		std::string get_summary_of_work_so_far() const ;
		std::size_t get_number_of_tasks_completed() const ;
		
	private:
		enum ThreadState { e_Idle = 0, e_Processing = 1, e_Terminating = 2 } ;
		std::vector< char > m_states ;
		std::vector< Task* > m_tasks ;
		std::vector< std::size_t > m_tasks_completed ;

		boost::thread_group	m_thread_group ;
	
		void worker_thread_loop( std::size_t thread_index ) ;

		// forbid copying, assignment.
		MultiThreadedWorker( MultiThreadedWorker const& other ) ;
		MultiThreadedWorker& operator=( MultiThreadedWorker const& other ) ;
	} ;

}

#endif
