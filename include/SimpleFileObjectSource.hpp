#ifndef __GTOOL_SimpleFileObjectSOURCE_HPP__
#define __GTOOL_SimpleFileObjectSOURCE_HPP__

#include "ObjectSource.hpp"
#include "FileUtil.hpp"
#include "Whitespace.hpp"

// Read objects from a file making no attempt to chunk reads
template< typename Object >
struct SimpleFileObjectSource: public ObjectSource< Object >
{
public:

	SimpleFileObjectSource( INPUT_FILE_PTR stream_ptr ) 
	: m_stream_ptr( stream_ptr )
	{}

	~SimpleFileObjectSource() {} ;

	void read( Object& object ) {
		assert( !check_if_empty() ) ;
		Whitespace whitespace ;
		(*m_stream_ptr) >> whitespace >> object >> whitespace ;
	}
	bool check_if_empty() {
		return !m_stream_ptr.get() || !m_stream_ptr->good() ;
	}
	operator bool() {
		return (*m_stream_ptr) ;
	}

protected:

	INPUT_FILE_PTR const& stream_ptr() { return m_stream_ptr ; }

private:

	INPUT_FILE_PTR m_stream_ptr ;
} ;

#endif
