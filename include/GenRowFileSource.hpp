#ifndef __GTOOL_GENROWFILESOURCE_HPP
#define __GTOOL_GENROWFILESOURCE_HPP

#include <string>
#include <vector>
#include <stddef.h>
#include "GenRow.hpp"
#include "SimpleFileObjectSource.hpp"
#include "ChainingFileObjectSource.hpp"
#include "FileUtil.hpp"
#include "genbin.hpp"

// Open a GEN input file, returning the appropriate type of Source object.
std::auto_ptr< ObjectSource< GenRow > > get_genrow_source_from_file( std::string filename ) ;
std::auto_ptr< ObjectSource< GenRow > > get_genrow_source_from_files( std::vector< std::string > filenames ) ;

typedef SimpleFileObjectSource< GenRow > SimpleTextFileGenRowSource ;

struct SimpleBinaryFileGenRowSource: public SimpleFileObjectSource< GenRow >
{
	typedef SimpleFileObjectSource< GenRow > base_t ;
	
	SimpleBinaryFileGenRowSource( INPUT_FILE_PTR a_stream_ptr )
		: base_t( a_stream_ptr )
	{
		genbin::read_offset( *stream_ptr(), &m_offset ) ;
		stream_ptr()->ignore( m_offset ) ;
		if( !(*stream_ptr())) {
			throw BadFileFormatException( "SimpleBinaryFileGenRowSource: unable to read (or to skip) offset - this can't be a valid binary gen file." ) ;
		}
	}

	SimpleBinaryFileGenRowSource& read( GenRow & row ) {
		row.read_from_binary_stream( *stream_ptr() ) ;
		return *this ;
	}
	
	private:
		uint32_t m_offset ;
} ;

struct ChainingFileGenRowSource: public ChainingFileObjectSource< GenRow >
{
	ChainingFileGenRowSource( std::vector< std::string > filenames ) ;
} ;

#endif
