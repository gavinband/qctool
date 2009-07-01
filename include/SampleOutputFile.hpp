#ifndef __GTOOL_SAMPLE_OUTPUT_FILE_HPP__
#define __GTOOL_SAMPLE_OUTPUT_FILE_HPP__

#include <vector>
#include <string>
#include "ObjectSink.hpp"
#include "SampleRow.hpp"
#include "GToolException.hpp"
#include "string_utils.hpp"

struct SampleOutputFileException: public GToolException
{
	SampleOutputFileException( std::string const& msg) : GToolException( msg ) {}
} ;

// A SampleRow sink representing a sample file, templated on
// the StreamSinkType so that different file sources can be used.
// The difference between this and the underlying StreamSinkType object
// is that, on the write of the first SampleRow, this writes the column
// header and column type lines at the beginning before acting like the
// StreamSinkType object.
// See http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format.html
// for a description of the sample file format.
template< typename StreamSinkType >
struct SampleOutputFile: public ObjectSink< SampleRow >
{
	SampleOutputFile( OUTPUT_FILE_PTR stream_ptr )
	: 	m_stream_ptr( stream_ptr ),
		m_have_written_header_and_types( false )
	{
		assert( m_stream_ptr.get() != 0 ) ;
	}

	SampleOutputFile& write( SampleRow const& row ) {
		if( !m_have_written_header_and_types ) {
			write_header_and_types( row.column_headings(), row.column_types() ) ;
			if( !*m_stream_ptr ) {
				throw SampleOutputFileException( "SampleOutputFile: error encountered writing header and type rows." ) ;
			}
			m_row_sink_ptr.reset( new StreamSinkType( m_stream_ptr )) ;
		}
		m_row_sink_ptr->write( row ) ;
		return *this ;
	}

	operator bool() const {
		if( m_have_written_header_and_types ) {
			return (*m_stream_ptr) ;
		}
		else {
			return (*m_row_sink_ptr);
		}
	}

private:

	void write_header_and_types( std::vector< std::string > const& column_headings, std::vector< char > column_types ) {
		for( std::vector< std::string >::const_iterator i = column_headings.begin(); i != column_headings.end(); ++i ) {
			if( i != column_headings.begin() )
				(*m_stream_ptr) << " " ;
			(*m_stream_ptr) << *i ;
		}
		(*m_stream_ptr) << "\n" ;
		for( std::vector< char >::const_iterator i = column_types.begin(); i != column_types.end(); ++i ) {
			if( i != column_types.begin() )
				(*m_stream_ptr) << " " ;
			(*m_stream_ptr) << *i ;			
		}
		(*m_stream_ptr) << "\n" ;
		m_have_written_header_and_types = true ;
	}

	OUTPUT_FILE_PTR m_stream_ptr ;
	std::auto_ptr< StreamSinkType > m_row_sink_ptr ;
	bool m_have_written_header_and_types ;
} ;

#endif
