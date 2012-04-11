
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef ISTREAM_AGGREGATOR_HPP
#define ISTREAM_AGGREGATOR_HPP

#include <iostream>
#include <memory>

namespace statfile {
	struct IstreamAggregator
	{
	protected:
	
		std::istream& stream() { return *m_stream_ptr ; }
		void set_stream( std::auto_ptr< std::istream > stream_ptr ) {
			m_stream_ptr = stream_ptr ;
			turn_on_ios_exceptions() ;
		}

		void reset_stream_to_start() {
			reset_stream_to( 0u ) ;
		}
		
		void reset_stream_to( std::streamoff position ) {
			m_stream_ptr->clear() ;
			assert( *m_stream_ptr ) ;
			m_stream_ptr->seekg( position, std::ios::beg ) ;
		}
		
		void turn_off_ios_exceptions() const {
			m_stream_ptr->exceptions( std::ios::badbit ) ;
		}

		void turn_on_ios_exceptions() const {
			m_stream_ptr->exceptions( std::ios::failbit | std::ios::badbit ) ;
		}
		
	private:
	
		std::auto_ptr< std::istream > m_stream_ptr ;	
	} ;
}

#endif