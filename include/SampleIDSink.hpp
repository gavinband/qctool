#ifndef SAMPLEIDSINK_HPP
#define SAMPLEIDSINK_HPP

#include <vector>
#include <string>
#include "ObjectSink.hpp"
#include "SampleRow.hpp"
#include "GToolException.hpp"
#include "string_utils/string_utils.hpp"

struct SampleIDSinkException: public std::exception
{
	char const* what() const throw() { return "SampleIDSinkException" ; }
} ;

// A SampleRow sink which just outputs sample ids, one per line,
// to the given stream.
struct SampleIDSink: public ObjectSink< SampleRow >
{
	SampleIDSink( OUTPUT_FILE_PTR stream_ptr )
	: 	m_stream_ptr( stream_ptr )
	{
		assert( m_stream_ptr.get() != 0 ) ;
	}

	SampleIDSink& write( SampleRow const& row ) {
		(*m_stream_ptr) << row.ID1() << " " << row.ID2() << "\n" ;
		return *this ;
	}

	operator bool() const {
		return (*m_stream_ptr) ;
	}

private:

	OUTPUT_FILE_PTR m_stream_ptr ;
} ;

#endif
