
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <cassert>
#include "../../config.hpp"
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
#include "appcontext/FileUtil.hpp"

namespace appcontext {
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
		bio::file_source source( filename, open_mode ) ;
		if( !source.is_open() ) {
			throw FileNotOpenedError( filename ) ;
		}
		stream_ptr->push( source ) ;
	#else
		if( compression_flags != e_None ) {
			throw FileException( "File compression requested.  Please recompile with boost support.") ;		
		}
		std::auto_ptr< std::ifstream > stream_ptr( new std::ifstream( filename.c_str(), open_mode )) ;
		if( !stream_ptr->is_open() ) {
			throw FileNotOpenedError( filename ) ;
		}
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
		bio::file_sink sink( filename, open_mode ) ;
		if( !sink.is_open() ) {
			throw FileNotOpenedError( filename ) ;
		}
		stream_ptr->push( sink ) ;
	#else
		if( compression_flags != e_None ) {
			throw FileException( "File compression requested.  Please recompile with boost support.") ;		
		}
		std::auto_ptr< std::ofstream > stream_ptr( new std::ofstream( filename.c_str(), open_mode )) ;
		if( !stream_ptr->is_open() ) {
			throw FileNotOpenedError( filename ) ;
		}
	#endif

		return OUTPUT_FILE_PTR( stream_ptr.release() ) ;
	}

	//
	FileCompressionType determine_file_compression( std::string const& filename ) {
		if( filename.find( ".gz" ) == ( filename.size()-3 )) {	
			return e_Gzip ;
		}
		else if( filename.find( ".bz2" ) == ( filename.size()-4 )) {
			return e_Bzip2 ;
		}
		else {
			return e_None ;
		}
	}

	FileModeType determine_file_mode( std::string const& filename ) {
		if( filename.find( ".bgen" ) != std::string::npos ) {
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

	bool exists( std::string const& filename ) {
	#if HAVE_BOOST_FILESYSTEM
		return BFS::exists( filename ) ;
	#else
		assert(0) ;
	#endif
	}

	bool is_regular( std::string const& filename ) {
	#if HAVE_BOOST_FILESYSTEM
		return BFS::is_regular( filename ) ;
	#else
		assert(0) ;
	#endif
	}

	void copy_file( std::string const& filename1, std::string const& filename2 ) {
	#if HAVE_BOOST_FILESYSTEM
		return BFS::copy_file( filename1, filename2 ) ;
	#else
		assert(0) ;
	#endif
	}

	void rename( std::string const& filename1, std::string const& filename2 ) {
	#if HAVE_BOOST_FILESYSTEM
		return BFS::rename( filename1, filename2 ) ;
	#else
		assert(0) ;
	#endif
	}

	std::string remove_extension_if_present( std::string const& filename ) {
	#if HAVE_BOOST_FILESYSTEM
		BFS::path path( filename ) ;
		return (path.parent_path() / path.stem()).string() ;
	#else
		assert(0) ;
	#endif
	}

	std::string remove_all_extensions( std::string const& filename ) {
	#if HAVE_BOOST_FILESYSTEM
		BFS::path path( filename ) ;
		while( path.extension() != "" ) {
			path = path.parent_path() / path.stem() ;
		}
		return path.filename() ;
	#else
		assert(0) ;
	#endif
	}
}
