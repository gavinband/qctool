#include <iostream>
#include <string>
#include <cassert>
#include "FileUtil.hpp"
#include "../config.hpp"
#ifdef HAVE_BOOST_IOSTREAMS
	#include <boost/iostreams/filtering_stream.hpp>
	#include <boost/iostreams/filter/gzip.hpp>
	#include <boost/iostreams/device/file.hpp>
	namespace bio = boost::iostreams ;
#else
	#include <fstream>
#endif

#if HAVE_BOOST_FILESYSTEM
	#include <boost/filesystem.hpp>
	namespace BFS = boost::filesystem ;
#endif

// Return a stream representing a given input file, optionally with gzip decompression
INPUT_FILE_PTR
open_file_for_input( std::string const& filename, int mode_flags ) {
	std::ios_base::openmode open_mode = std::ios_base::in ;

	if(( mode_flags & e_FileModeMask ) == e_BinaryMode ) {
		open_mode |= std::ios_base::binary ;
	}

	int compression_flags = mode_flags & e_FileCompressionMask ;

#ifdef HAVE_BOOST_IOSTREAMS
	std::auto_ptr< bio::filtering_istream > stream_ptr( new bio::filtering_istream ) ;
	switch( compression_flags ) {
		case e_Gzip :
			stream_ptr->push(bio::gzip_decompressor());
			break ;
		case e_Bzip2 :
			assert( 0 ) ;			// Bzip2 not yet supported.
			break ;
		default:
			break ;
	}
	stream_ptr->push( bio::file_source( filename, open_mode )) ;
#else
	if( compression_flags != e_None ) {
		throw FileException( "File compression requested.  Please recompile with boost support.") ;		
	}
	std::auto_ptr< std::ifstream > stream_ptr( new std::ifstream( filename.c_str(), open_mode )) ;
#endif
	
	return INPUT_FILE_PTR( stream_ptr.release() ) ;
}

// Return a stream representing a given output file, optionally with gzip compression.
OUTPUT_FILE_PTR
open_file_for_output( std::string const& filename, int mode_flags ) {
	std::ios_base::openmode open_mode = std::ios_base::out ;
	if(( mode_flags & e_FileModeMask ) == e_BinaryMode ) {
		open_mode |= std::ios_base::binary ;
	}

	int compression_flags = mode_flags & e_FileCompressionMask ;

#ifdef HAVE_BOOST_IOSTREAMS
	std::auto_ptr< bio::filtering_ostream > stream_ptr( new bio::filtering_ostream ) ;
	switch( compression_flags ) {
		case e_Gzip:
			stream_ptr->push(bio::gzip_compressor());
			break ;
		case e_Bzip2:
			assert( 0 ) ; // not supported yet.
			break ;
		case e_None:
			break ;
		default:
			assert(0) ; 
			break ;
	}
	stream_ptr->push( bio::file_sink( filename, open_mode )) ;
#else
	if( compression_flags != e_None ) {
		throw FileException( "File compression requested.  Please recompile with boost support.") ;		
	}
	std::auto_ptr< std::ofstream > stream_ptr( new std::ofstream( filename.c_str(), open_mode )) ;
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

FileModeType determine_file_mode( std::string const& filename ) {
	if( filename.find( ".bin" ) != std::string::npos ) {
		return e_BinaryMode ;
	}
	else {
		return e_TextMode ;
	}
}


// Return a stream representing a given input file, attempting to auto-detect any compression used.
INPUT_FILE_PTR
open_file_for_input( std::string const& filename ) {
	FileCompressionType file_compression = determine_file_compression( filename ) ;
	FileModeType file_mode = determine_file_mode( filename ) ;
	return open_file_for_input( filename, file_compression | file_mode ) ;
}


// Return a stream representing a given input file, attempting to auto-detect the compression to
// use from the filename.
OUTPUT_FILE_PTR
open_file_for_output( std::string const& filename ) {
	FileCompressionType file_compression = determine_file_compression( filename ) ;
	FileModeType file_mode = determine_file_mode( filename ) ;
	return open_file_for_output( filename, file_compression | file_mode ) ;
}

std::vector< std::string > find_files_matching_path_with_wildcard( std::string filename_with_wildcard ) {
	std::vector< std::string > result ;
	std::size_t wildcard_pos = filename_with_wildcard.find( '*' ) ;
#if HAVE_BOOST_FILESYSTEM
	BFS::path dir = BFS::path( filename_with_wildcard ).parent_path() ;
  	filename_with_wildcard = BFS::path( filename_with_wildcard ).filename() ;
	if( dir.empty() ) {
		dir = "." ;
	}
	BFS::directory_iterator dir_i( dir ), end_i ;

	std::size_t wildcard_pos = filename_with_wildcard.find( '*' ) ;
	std::string first_bit = filename_with_wildcard.substr( 0, wildcard_pos ) ;
	std::string second_bit = "" ;
	if( wildcard_pos != std::string::npos ) {
		second_bit = filename_with_wildcard.substr( wildcard_pos + 1, filename_with_wildcard.size()) ;
	}
	
	for( ; dir_i != end_i; ++dir_i ) {
		if( BFS::is_regular_file( *dir_i )) {
			std::string filename = dir_i->filename();
			if( filename.size() >= ( first_bit.size() + second_bit.size()) ) {
				if( first_bit == filename.substr( 0, first_bit.size()) ) {
					if(( filename.size() == first_bit.size())
					|| (( wildcard_pos != std::string::npos )
					&& (second_bit == filename.substr( filename.size() - second_bit.size(), second_bit.size())))) {
						result.push_back(( dir / filename ).string() ) ;
					}
				}
			}
		}
	}
#else
	// Oh dear.  Assume the filename is valid as is.
	result.push_back( filename_with_wildcard ) ;
#endif
	return result ;
}

bool exists( std::string const& filename ) {
#if HAVE_BOOST_FILESYSTEM
	return BFS::exists( filename ) ;
#endif
}

bool is_regular( std::string const& filename ) {
#if HAVE_BOOST_FILESYSTEM
	return BFS::is_regular( filename ) ;
#endif
}
