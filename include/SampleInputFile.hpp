#ifndef __GTOOL_SAMPLE_INPUT_FILE_HPP__
#define __GTOOL_SAMPLE_INPUT_FILE_HPP__

#include <vector>
#include <string>
#include "ObjectSource.hpp"
#include "SampleRow.hpp"
#include "GToolException.hpp"
#include "string_utils.hpp"

struct SampleInputFileException: public GToolException
{
	SampleInputFileException( std::string const& msg) : GToolException( msg ) {}
} ;

// A SampleRow source representing a sample file, templated on
// the FileSourceType so that different file sources can be used.
// The difference between this and a plain FileSourceType object
// is that this reads the column header and column type lines at the beginning
// of the sample file before acting like the FileSourceType object.
// See http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format.html
// for a description of the sample file format.
template< typename FileSourceType >
struct SampleInputFile: public ObjectSource< SampleRow >
{
	SampleInputFile( INPUT_FILE_PTR stream_ptr )
	{
		assert( stream_ptr.get() != 0 ) ;
		// Read the column headings and column row.
		std::string line ;
		std::getline( *stream_ptr, line ) ;
		m_column_headers = split_discarding_empty_entries( line, " " ) ;
		std::getline( *stream_ptr, line ) ;
		std::vector<std::string> column_type_vector = split_discarding_empty_entries( line, " " ) ;
		m_column_types.resize( column_type_vector.size() ) ;
		for( std::size_t i = 0; i < column_type_vector.size(); ++i ) {
			m_column_types[i] = column_type_vector[i][0] ;
		}
		if( !stream_ptr->good()) {
			throw SampleInputFileException( "A problem occured during or after reading the column headings and types." ) ;
		}
		
		// Construct the source, handing over our stream pointer.
		m_row_source_ptr.reset( new FileSourceType( stream_ptr )) ;
	}

	// required methods from ObjectSource, these forward
	// to the SourceType member.
	void read( SampleRow& row ) {
		row.reset( m_column_headers, m_column_types ) ;
		m_row_source_ptr->read( row ) ;
	}
	bool check_if_empty() { return !m_row_source_ptr.get() || m_row_source_ptr->check_if_empty() ; }

	operator bool() { return *m_row_source_ptr ; }
private:

	std::auto_ptr< FileSourceType > m_row_source_ptr ;
	std::vector< std::string > m_column_headers ;
	std::vector< char > m_column_types ;
} ;

#endif
