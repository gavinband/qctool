#ifndef TIMER_HPP
#define TIMER_HPP

#include "../config.hpp"
#if HAVE_BOOST_TIMER
	#include <boost/timer.hpp>
	typedef boost::timer Timer ;
#else
	struct Timer
	{
		double elapsed_time() const { return 0.0 ; }
		void restart() {} ;
	} ;
#endif


#endif
