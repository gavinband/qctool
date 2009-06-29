#ifndef __GTOOL_GENROWFILESINK_HPP
#define __GTOOL_GENROWFILESINK_HPP

#include "GenRow.hpp"
#include "SimpleFileObjectSink.hpp"
#include "FileUtil.hpp"
typedef SimpleFileObjectSink< GenRow > SimpleGenRowTextFileSink ;

struct SimpleGenRowBinaryFileSink: public SimpleFileObjectSink< GenRow >
{
	typedef SimpleFileObjectSink< GenRow > base_t ;

	SimpleGenRowBinaryFileSink( OUTPUT_FILE_PTR stream_ptr )
		: base_t( stream_ptr )
	{}

	void write( GenRow const& row ) {
		assert( !check_if_full() ) ;
		row.write_to_binary_stream( *stream_ptr() ) ;
	}
} ;

#endif
