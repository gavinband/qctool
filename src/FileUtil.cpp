#include <iostream>
#include <string>
#include "FileUtil.hpp"
#include "../config.hpp"
#ifdef HAVE_BOOST_IOSTREAMS
	#include <boost/iostreams/filtering_stream.hpp>
	#include <boost/iostreams/filter/gzip.hpp>
	#include <boost/iostreams/device/file.hpp>
	namespace bio = boost::iostreams ;
#endif

// Return a stream representing a given input file, optionally with gzip decompression
INPUT_FILE_PTR
open_file_for_input( std::string const& filename, FileCompressionType file_compression_type ) {
#ifdef HAVE_BOOST_IOSTREAMS
	std::auto_ptr< bio::filtering_istream > stream_ptr( new bio::filtering_istream ) ;
	switch( file_compression_type == e_Gzip ) {
		case e_Gzip :
			stream_ptr->push(bio::gzip_decompressor());
			break ;
		case e_Bzip2 :
			assert( 0 ) ;			// Bzip2 not yet supported.
			break ;
		default:
			break ;
	}
	stream_ptr->push( bio::file_source( filename )) ;
#else
	if( file_compression_type != e_None ) {
		throw FileException( "File compression requested.  Please recompile with boost support.") ;		
	}
	std::auto_ptr< std::ifstream > stream_ptr( new std::ifstream( filename )) ;
#endif
	
	return INPUT_FILE_PTR( stream_ptr.release() ) ;
}

// Return a stream representing a given output file, optionally with gzip compression.
OUTPUT_FILE_PTR
open_file_for_output( std::string const& filename, FileCompressionType file_compression_type ) {
#ifdef HAVE_BOOST_IOSTREAMS
	std::auto_ptr< bio::filtering_ostream > stream_ptr( new bio::filtering_ostream ) ;
	switch( file_compression_type ) {
		case e_Gzip:
			stream_ptr->push(bio::gzip_compressor());
			break ;
		case e_Bzip2:
			assert( 0 ) ; // not supported yet.
			break ;
		default:
			break ;
	}
	stream_ptr->push( bio::file_sink( filename )) ;
#else
	if( file_compression_type != e_None ) {
		throw FileException( "File compression requested.  Please recompile with boost support.") ;		
	}
	std::auto_ptr< std::ifstream > stream_ptr( new std::ifstream( filename )) ;
#endif
	return OUTPUT_FILE_PTR( stream_ptr.release() ) ;
}

//
FileCompressionType determine_file_compression( std::string const& filename ) {
	if( filename.find( ".gz" ) != std::string::npos ) {
		return e_Gzip ;
	}
	else if( filename.find( ".bz2" ) != std::string::npos ) {
		return e_Bzip2 ;
	}
	else {
		return e_None ;
	}
}


// Return a stream representing a given input file, attempting to auto-detect any compression used.
INPUT_FILE_PTR
open_file_for_input( std::string const& filename ) {
	FileCompressionType file_compression = determine_file_compression( filename ) ;
	return open_file_for_input( filename, file_compression ) ;
}


// Return a stream representing a given input file, attempting to auto-detect the compression to
// use from the filename.
OUTPUT_FILE_PTR
open_file_for_output( std::string const& filename ) {
	return open_file_for_output( filename, determine_file_compression( filename )) ;
}


