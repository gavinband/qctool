#ifndef TIMER_HPP
#define TIMER_HPP

#include <string>
#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "../../config.hpp"


#if HAVE_MACH_TIME
#include <mach/mach_time.h>
namespace appcontext {
	struct Timer
	{
		Timer() {
			mach_timebase_info_data_t timebase_info ;
			kern_return_t err = mach_timebase_info( &timebase_info ) ;
			assert( err == 0 ) ;
			m_factor = (double( timebase_info.numer ) / double( timebase_info.denom )) / 1000000000.0 ;
			restart() ;
		}

		double elapsed() const {
			return double(mach_absolute_time() - m_start_time) * m_factor ;
		}

		void restart() {
			m_start_time = mach_absolute_time() ;
		} ;

		std::string display() const {
			std::ostringstream os ;
			os << std::fixed << std::setprecision(1) << elapsed() << "s" ;
			return os.str() ;
		}

	private:
		
		uint64_t m_start_time ;
		double m_factor ;
		mach_timebase_info_data_t m_mach_timebase_info ;
	} ;
}
#elif HAVE_BOOST_TIMER
#include <boost/timer.hpp>
namespace appcontext {
	struct Timer
	{
		double elapsed() const {
			return m_boost_timer.elapsed() ;
		}
		
		void restart() {
			m_boost_timer.restart() ;
		}
		
		std::string display() const {
			std::ostringstream os ;
			os << std::fixed << std::setprecision(1) << elapsed() << "s" ;
			return os.str() ;
		}

	private:
		
		boost::timer m_boost_timer ;
	} ;
}
#elif HAVE_GETTIMEOFDAY
#include <sys/time.h>
namespace appcontext {
	struct Timer
	{
		Timer()
		{
			gettimeofday( &m_start_time, NULL ) ;
		}
		
		double elapsed() const {
			timeval current_time ;
			gettimeofday( &current_time, NULL ) ;
			return(
				double( current_time.tv_sec - m_start_time.tv_sec )
				+ (double( current_time.tv_usec ) - double( m_start_time.tv_usec )) / 100000.0
			) ;
		}

		void restart() {
			gettimeofday( &m_start_time, NULL ) ;
		}

		std::string display() const {
			std::ostringstream os ;
			os << std::fixed << std::setprecision(1) << elapsed() << "s" ;
			return os.str() ;
		}

	private:
		timeval m_start_time ;
	} ;
}
#else
namespace appcontext {
	struct Timer
	{
		double elapsed() const { return 0.0 ; }
		void restart() {} ;
		std::string display() const { return "?s "; }
	} ;
}
#endif

#endif
