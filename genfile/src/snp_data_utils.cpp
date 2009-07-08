#include <iostream>
#include <string>
#include <cstdio>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include "snp_data_utils.hpp"

namespace genfile {
	bool filename_indicates_bgen_format( std::string const& filename ) {
		return ( filename.find( ".bgen") != std::string::npos ) ;
	}

	CompressionType get_compression_type_indicated_by_filename( std::string const& filename ) {
		if( filename.find( ".gz") != std::string::npos ) {
			return e_GzipCompression ;
		}
		else {
			return e_NoCompression ;
		}
	}

	std::auto_ptr< std::istream > open_text_file_for_input( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_istream > gen_file_ptr( new boost::iostreams::filtering_istream ) ;
	    if (compression_type == e_GzipCompression) gen_file_ptr->push(boost::iostreams::gzip_decompressor());
		gen_file_ptr->push(boost::iostreams::file_source(filename.c_str())); 
		return std::auto_ptr< std::istream >( gen_file_ptr ) ;
  	}

	std::auto_ptr< std::ostream > open_text_file_for_output( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_ostream > gen_file_ptr( new boost::iostreams::filtering_ostream ) ;
	    if (compression_type == e_GzipCompression) gen_file_ptr->push(boost::iostreams::gzip_compressor());
		gen_file_ptr->push(boost::iostreams::file_sink(filename.c_str())); 
		return std::auto_ptr< std::ostream >( gen_file_ptr ) ;
  	}

	std::auto_ptr< std::istream > open_binary_file_for_input( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_istream > gen_file_ptr( new boost::iostreams::filtering_istream ) ;
	    if (compression_type == e_GzipCompression) gen_file_ptr->push(boost::iostreams::gzip_decompressor()) ;
		gen_file_ptr->push( boost::iostreams::file_source( filename.c_str(), std::ios::binary )) ;
		return std::auto_ptr< std::istream >( gen_file_ptr ) ;
  	}

	std::auto_ptr< std::ostream > open_binary_file_for_output( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_ostream > gen_file_ptr( new boost::iostreams::filtering_ostream ) ;
	    if (compression_type == e_GzipCompression) gen_file_ptr->push(boost::iostreams::gzip_compressor()) ;
		gen_file_ptr->push( boost::iostreams::file_sink( filename.c_str(), std::ios::binary )) ;
		return std::auto_ptr< std::ostream >( gen_file_ptr ) ;
  	}

	std::string create_temporary_filename() {
		return std::tmpnam( 0 ) ;
	}

	Ignorer ignore() { return Ignorer() ; }	
}
