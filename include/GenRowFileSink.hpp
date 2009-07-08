#ifndef __GTOOL_GENROWFILESINK_HPP
#define __GTOOL_GENROWFILESINK_HPP

#include "GenRow.hpp"
#include "SimpleFileObjectSink.hpp"
#include "FileUtil.hpp"
#include "SNPDataSink.hpp"

typedef SimpleFileObjectSink< GenRow > SimpleGenRowTextFileSink ;

struct SimpleGenRowBinaryFileSink: public SimpleFileObjectSink< GenRow >
{
	typedef SimpleFileObjectSink< GenRow > base_t ;

	SimpleGenRowBinaryFileSink( OUTPUT_FILE_PTR a_stream_ptr )
		: base_t( a_stream_ptr )
	{
		genfile::bgen::uint32_t offset = 0 ;
		genfile::bgen::write_offset( *stream_ptr(), offset ) ;
	}

	SimpleGenRowBinaryFileSink& write( GenRow const& row ) {
		row.write_to_binary_stream( *stream_ptr() ) ;
		return *this ;
	}
} ;

struct SNPDataSinkGenRowSink: public ObjectSink< GenRow >
{
public:
	SNPDataSinkGenRowSink( std::auto_ptr< genfile::SNPDataSink > snp_data_sink ) ;
	operator bool() const { return *m_snp_data_sink ; }
	SNPDataSinkGenRowSink& write( GenRow const& row ) ;
	
protected:

		std::auto_ptr< genfile::SNPDataSink > m_snp_data_sink ;
	
} ;

#endif
