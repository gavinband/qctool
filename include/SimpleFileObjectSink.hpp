#ifndef __GTOOL_SimpleFileObjectSink_HPP__
#define __GTOOL_SimpleFileObjectSink_HPP__

#include "ObjectSink.hpp"
#include "FileUtil.hpp"

// Read gen rows from a gen file making no attempt to chunk reads
template< typename Object >
struct SimpleFileObjectSink: public ObjectSink< Object >
{
public:

	SimpleFileObjectSink( OUTPUT_FILE_PTR stream_ptr )
	: m_stream_ptr( stream_ptr )
	{}

	void write( Object const& object ) {
		assert( !check_if_full() ) ;
		(*m_stream_ptr) << object ;
	}

	bool check_if_full() {
		return !m_stream_ptr.get() || !m_stream_ptr->good() ;
	}

protected:

	OUTPUT_FILE_PTR const& stream_ptr() { return m_stream_ptr ; }

private:

	OUTPUT_FILE_PTR m_stream_ptr ;
} ;

#endif