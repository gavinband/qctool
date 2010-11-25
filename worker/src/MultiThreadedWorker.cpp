#include <sstream>
#include <iomanip>
#include <numeric>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include "worker/MultiThreadedWorker.hpp"

namespace worker
{
	MultiThreadedWorker::MultiThreadedWorker( std::size_t number_of_threads )
		: m_states( number_of_threads, e_Idle ),
		  m_tasks( number_of_threads, 0 ),
		  m_tasks_completed( number_of_threads, 0 )
	{
		for( std::size_t i = 0; i < m_states.size(); ++i ) {
			m_thread_group.create_thread(
				boost::bind(
					&MultiThreadedWorker::worker_thread_loop,
					this,
					i
				)
			) ;
		}
	}


	MultiThreadedWorker::~MultiThreadedWorker() {
		for( std::size_t i = 0; i < m_states.size(); ++i ) {
			while( m_states[i] == e_Processing ) {
				boost::this_thread::yield() ;
			}
			m_states[i] = e_Terminating ;
		}
		m_thread_group.join_all() ;		
	}

	bool MultiThreadedWorker::ask_to_perform_task( Task& task ) {
		for( std::size_t i = 0; i < m_states.size(); ++i ) {
			if( m_states[i] == e_Idle ) {
				m_tasks[i] = &task ;
				m_states[i] = e_Processing ;
				return true ;
			}
		}
		return false ;
	}

	void MultiThreadedWorker::tell_to_perform_task( Task& task ) {
		while( !ask_to_perform_task( task )) {
			boost::this_thread::sleep( boost::posix_time::milliseconds( 10 ) ) ;
		}
	}

	std::size_t MultiThreadedWorker::get_number_of_worker_threads() const {
		return m_thread_group.size() ;
	}

	std::size_t MultiThreadedWorker::get_number_of_tasks_performed_by_thread( std::size_t index ) const {
		assert( index < m_tasks_completed.size() ) ;
		return m_tasks_completed[ index ] ;
	}

	std::string MultiThreadedWorker::get_summary_of_work_so_far() const {
		std::ostringstream oStream ;
		oStream << "Tasks carried out by " << get_number_of_worker_threads() << " worker threads:\n" ;
		for( std::size_t i = 0; i < get_number_of_worker_threads(); ++i ) {
			oStream << "thread " << std::setw(3) << i << ": " << get_number_of_tasks_performed_by_thread( i ) << "\n" ;
		}
		return oStream.str() ;
	}

	std::size_t MultiThreadedWorker::get_number_of_tasks_completed() const {
		return std::accumulate( m_tasks_completed.begin(), m_tasks_completed.end(), 0 ) ;
	}

	void MultiThreadedWorker::worker_thread_loop( std::size_t thread_index ) {
		char& state = m_states[ thread_index ] ;
		while( 1 ) {
			switch( state ) {
				case e_Idle:
					boost::this_thread::yield() ;
					break ;
				case e_Terminating:
					return ;
					break ;
				case e_Processing:
					perform_task( *m_tasks[ thread_index ] ) ;
					++m_tasks_completed[ thread_index ] ;
					state = e_Idle ;
					break ;
			}
		}
	}
}
