#ifndef __GTOOL_GENROWFILESOURCE_HPP
#define __GTOOL_GENROWFILESOURCE_HPP

#include <string>
#include <vector>
#include "GenRow.hpp"
#include "SimpleFileObjectSource.hpp"
#include "ChainingFileObjectSource.hpp"
#include "FileUtil.hpp"

// Open a GEN input file, returning the appropriate type of Source object.
std::auto_ptr< ObjectSource< GenRow > > get_genrow_source_from_file( std::string filename ) ;
std::auto_ptr< ObjectSource< GenRow > > get_genrow_source_from_files( std::vector< std::string > filenames ) ;

typedef SimpleFileObjectSource< GenRow > SimpleTextFileGenRowSource ;

struct SimpleBinaryFileGenRowSource: public SimpleFileObjectSource< GenRow >
{
	typedef SimpleFileObjectSource< GenRow > base_t ;
	
	SimpleBinaryFileGenRowSource( INPUT_FILE_PTR stream_ptr )
		: base_t( stream_ptr )
	{}

	void read( GenRow & row ) {
		assert( !check_if_empty() ) ;
		row.read_from_binary_stream( *stream_ptr() ) ;
	}
} ;

struct ChainingFileGenRowSource: public ChainingFileObjectSource< GenRow >
{
	ChainingFileGenRowSource( std::vector< std::string > filenames ) ;
} ;

#endif
