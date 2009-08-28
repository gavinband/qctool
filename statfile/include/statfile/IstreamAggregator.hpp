#ifndef ISTREAM_AGGREGATOR_HPP
#define ISTREAM_AGGREGATOR_HPP

#include <iostream>
#include <memory>

namespace statfile {
	struct IstreamAggregator
	{
	public:
	
		operator void*() const { return *m_stream_ptr ; }

	protected:
	
		std::istream& stream() { return *m_stream_ptr ; }
		void set_stream( std::auto_ptr< std::istream > stream_ptr ) { m_stream_ptr = stream_ptr ; }

		void reset_stream_to_start() {
			m_stream_ptr->clear() ;
			assert( *m_stream_ptr ) ;
			m_stream_ptr->seekg(0) ;
		}
		
	private:
	
		std::auto_ptr< std::istream > m_stream_ptr ;	
	} ;
}

#endif