
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef OSTREAM_STATSINK_HPP
#define OSTREAM_STATSINK_HPP

#include <iostream>
#include <memory>

namespace statfile {
	struct OstreamAggregator
	{
	public:
	
		operator bool() const { return (*m_stream_ptr) ; }

	protected:
	
		std::ostream& stream() { return *m_stream_ptr ; }
		void set_stream( std::auto_ptr< std::ostream > stream_ptr ) { m_stream_ptr = stream_ptr ; }

	private:
	
		std::auto_ptr< std::ostream > m_stream_ptr ;	
	} ;
}

#endif