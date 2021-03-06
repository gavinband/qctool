
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef WORKER_FUNCTION_TASK_HPP
#define WORKER_FUNCTION_TASK_HPP

#include <boost/function.hpp>
#include "worker/Task.hpp"

namespace worker {
	class FunctionTask: public Task
	{
	public:
		template< typename Function >
		FunctionTask( Function const& function ):
			m_function( function )
		{
			assert( m_function ) ;
		}
		
		void operator()() {
			m_function() ;
		}
		
	private:
		boost::function< void() > m_function ;
		FunctionTask( FunctionTask const& other ) ; //forbid copying.
	} ;
}

#endif
