#include <boost/thread.hpp>
#include "test_case.hpp"
#include "worker/Worker.hpp"
#include "worker/SynchronousWorker.hpp"
#include "worker/MultiThreadedWorker.hpp"
#include "worker/QueuedMultiThreadedWorker.hpp"

struct Multiply: public worker::Task
{
	Multiply(): m_value1(0), m_value2(0) {}

	Multiply( double& value1, double const& value2  )
		: m_value1( &value1 ),
		  m_value2( &value2 )
	{}

	Multiply( Multiply const& other ) 
		: m_value1( other.m_value1 ),
		  m_value2( other.m_value2 )
	{}

	Multiply& operator=( Multiply const& other ) {
		m_value1 = other.m_value1 ;
		m_value2 = other.m_value2 ;
		return *this ;
	}
	
	void operator()() {
		(*m_value1) *= (*m_value2) ;
	}
	
private:
	
	double* m_value1 ;
	double const* m_value2 ;
} ;

struct TestFailedException: public std::exception
{
	TestFailedException( worker::Task* task, double expected, double got )
		: m_task_address( task ),
		  m_expected( expected ),
		  m_got( got )
	{}
	
	~TestFailedException() throw() {}
	
	char const* what() const throw() { return "TestFailedException" ; }

	worker::Task const* get_task_address() const { return m_task_address ; }
	double get_expected_value() const { return m_expected ; }
	double get_actual_value() const { return m_got ; }
	
private:
	
	worker::Task const* m_task_address ;
	double m_expected ;
	double m_got ;
} ;

void test_worker( worker::Worker& worker ) {
	std::vector< double > numbers( 500 ) ;
	std::vector< Multiply > tasks ;
	tasks.reserve( numbers.size() ) ;
	for( std::size_t i = 0; i < numbers.size(); ++i ) {
		numbers[i] = double(i) ;
		tasks.push_back( Multiply( numbers[i], 2.0 )) ;
	}

	// Ask the worker to do all the tasks
	for( std::size_t i = 0; i < tasks.size(); ++i ) {
		while( !worker.ask_to_perform_task( tasks[i] )) {
			boost::this_thread::yield() ;
		}
	}

	// Wait until the worker has finished
	while( worker.get_number_of_tasks_completed() != numbers.size() ) {
		boost::this_thread::yield() ;
	}

	// Check all of the multiplications were performed
	for( std::size_t i = 0; i < numbers.size(); ++i ) {
		double expected = double(i) * 2.0 ;
		if( numbers[i] != expected ) {
			throw TestFailedException( &tasks[i], expected, numbers[i] ) ;
		}
	}
}

void test_worker2( worker::Worker& worker ) {
	std::vector< double > numbers( 500 ) ;
	std::vector< Multiply > tasks ;
	tasks.reserve( numbers.size() ) ;
	for( std::size_t i = 0; i < numbers.size(); ++i ) {
		numbers[i] = double(i) ;
		tasks.push_back( Multiply( numbers[i], 2.0 )) ;
	}

	// Ask the worker to do all the tasks
	for( std::size_t i = 0; i < tasks.size(); ++i ) {
		worker.tell_to_perform_task( tasks[i] ) ;
	}

	// Wait until the worker has finished
	while( worker.get_number_of_tasks_completed() != numbers.size() ) {
		boost::this_thread::yield() ;
	}

	// Check all of the multiplications were performed
	for( std::size_t i = 0; i < numbers.size(); ++i ) {
		double expected = double(i) * 2.0 ;
		if( numbers[i] != expected ) {
			throw TestFailedException( &tasks[i], expected, numbers[i] ) ;
		}
	}
}


AUTO_TEST_CASE( test_synchronous_worker ) {
	std::cerr << "Testing synchronous worker...\n" ;
	try {
		worker::SynchronousWorker worker ;
		test_worker( worker ) ;
		test_worker2( worker ) ;
	}
	catch( TestFailedException const& e ) {
		std::cerr << "!! Error: task " << e.get_task_address() << ": expected " << e.get_expected_value() << ", got: " << e.get_actual_value() << ".\n" ;
		TEST_ASSERT(0) ;
	}
}
AUTO_TEST_CASE( test_multi_threaded_worker ) {
	std::cerr << "Testing multi-threaded worker...\n" ;
	for( std::size_t number_of_threads = 1; number_of_threads < 100; ++number_of_threads ) {
		try {
			worker::MultiThreadedWorker worker( number_of_threads ) ;
			test_worker( worker ) ;
			test_worker2( worker ) ;
		}
		catch( TestFailedException const& e ) {
			std::cerr << "!! Error: task " << e.get_task_address() << ": expected " << e.get_expected_value() << ", got: " << e.get_actual_value() << ".\n" ;
			TEST_ASSERT(0) ;
		}
	}
}

AUTO_TEST_CASE( test_queued_multi_threaded_worker ) {
	std::cerr << "Testing queued multi-threaded worker...\n" ;
	for( std::size_t number_of_threads = 1; number_of_threads < 100; ++number_of_threads ) {
		std::cerr << "Number of threads = " << number_of_threads << ".\n" ;
		try {
			{
				worker::QueuedMultiThreadedWorker worker( number_of_threads ) ;
				test_worker( worker ) ;
				std::cerr << "Completed test, destructing worker...\n" ;
			}
			{
				worker::QueuedMultiThreadedWorker worker( number_of_threads ) ;
				test_worker2( worker ) ;
				std::cerr << "Completed test, destructing worker...\n" ;
			}
		}
		catch( TestFailedException const& e ) {
			std::cerr << "!! Error: task " << e.get_task_address() << ": expected " << e.get_expected_value() << ", got: " << e.get_actual_value() << ".\n" ;
			TEST_ASSERT(0) ;
		}
	}
}

AUTO_TEST_MAIN {
	test_synchronous_worker() ;
	test_multi_threaded_worker() ;
	test_queued_multi_threaded_worker() ;
}
