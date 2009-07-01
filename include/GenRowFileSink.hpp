#ifndef __GTOOL_GENROWFILESINK_HPP
#define __GTOOL_GENROWFILESINK_HPP

#include "GenRow.hpp"
#include "SimpleFileObjectSink.hpp"
#include "FileUtil.hpp"
#include "genbin.hpp"

typedef SimpleFileObjectSink< GenRow > SimpleGenRowTextFileSink ;

struct SimpleGenRowBinaryFileSink: public SimpleFileObjectSink< GenRow >
{
	typedef SimpleFileObjectSink< GenRow > base_t ;

	SimpleGenRowBinaryFileSink( OUTPUT_FILE_PTR stream_ptr )
		: base_t( stream_ptr )
	{
		genbin::write_offset( *stream_ptr, 0 ) ;
	}

	SimpleGenRowBinaryFileSink& write( GenRow const& row ) {
		row.write_to_binary_stream( *stream_ptr() ) ;
		return *this ;
	}
} ;

#endif
