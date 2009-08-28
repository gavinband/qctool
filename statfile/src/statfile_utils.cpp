#include <iostream>
#include <string>
#include <cstdio>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include "statfile/statfile_utils.hpp"

namespace statfile {
	IgnoreOne ignore() { return IgnoreOne() ; }
	IgnoreAll ignore_all() { return IgnoreAll() ; }
	EndRow end_row() { return EndRow() ; }

	ColumnSpec::ColumnSpec( std::string name, std::size_t repeats ) 
		: m_name( name ), m_number_of_repeats( repeats )
	{
		assert( m_number_of_repeats > 0 ) ;
	}

	std::vector< ColumnSpec > operator|( ColumnSpec const& left, ColumnSpec const& right ) {
		std::vector< ColumnSpec > result ;
		result.reserve(2) ;
		result.push_back( left ) ;
		result.push_back( right ) ;
		return result ;
	}

	std::vector< ColumnSpec > operator|( std::vector< ColumnSpec > const& left, ColumnSpec const& right ) {
		std::vector< ColumnSpec > result = left ;
		result.push_back( right ) ;
		return result ;
	}
	
	FileFormatType get_file_format_type_indicated_by_filename( std::string const& filename ) {
		if( filename.size() > 6 && filename.substr( filename.size() - 6, 6 ) == ".rdata" ) {
			return e_RFormat ;
		}
		return e_UnknownFormat ;
	}

	CompressionType get_compression_type_indicated_by_filename( std::string const& filename ) {
		if( filename.find( ".gz") != std::string::npos ) {
			return e_GzipCompression ;
		}
		else {
			return e_NoCompression ;
		}
	}

	std::string strip_file_format_extension_if_present( std::string const& filename ) {
		std::string recognised_extensions[1] = {
			".rdata"
		} ;

		for( std::size_t i = 0; i < 1u; ++i ) {
			if( filename.substr( filename.size() - recognised_extensions[i].size(), recognised_extensions[i].size() ) == recognised_extensions[i] ) {
				return filename.substr( 0, filename.size() - recognised_extensions[i].size() ) ;
			}
		}

		return filename ;
	}

	std::auto_ptr< std::istream > open_text_file_for_input( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_istream > gen_file_ptr( new boost::iostreams::filtering_istream ) ;
	    if (compression_type == e_GzipCompression) gen_file_ptr->push(boost::iostreams::gzip_decompressor());
		boost::iostreams::file_source file( filename.c_str() ) ;
		if( !file.is_open() ) {
			throw FileNotOpenedError() ;
		}
		gen_file_ptr->push( file ) ;
		return std::auto_ptr< std::istream >( gen_file_ptr ) ;
  	}

	std::auto_ptr< std::ostream > open_text_file_for_output( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_ostream > gen_file_ptr( new boost::iostreams::filtering_ostream ) ;
	    if (compression_type == e_GzipCompression) gen_file_ptr->push(boost::iostreams::gzip_compressor());
		boost::iostreams::file_sink file(filename.c_str()) ;
		if( !file.is_open() ) {
			throw FileNotOpenedError() ;
		}
		gen_file_ptr->push( file ); 
		return std::auto_ptr< std::ostream >( gen_file_ptr ) ;
  	}

	std::auto_ptr< std::istream > open_binary_file_for_input( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_istream > gen_file_ptr( new boost::iostreams::filtering_istream ) ;
	    if (compression_type == e_GzipCompression) gen_file_ptr->push( boost::iostreams::gzip_decompressor() ) ;
		boost::iostreams::file_source file( filename.c_str(), std::ios::binary ) ;
		if( !file.is_open() ) {
			throw FileNotOpenedError() ;
		}
		gen_file_ptr->push( file ) ;
		return std::auto_ptr< std::istream >( gen_file_ptr ) ;
  	}

	std::auto_ptr< std::ostream > open_binary_file_for_output( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_ostream > gen_file_ptr( new boost::iostreams::filtering_ostream ) ;
	    if (compression_type == e_GzipCompression) gen_file_ptr->push(boost::iostreams::gzip_compressor()) ;
		boost::iostreams::file_sink file(filename.c_str(), std::ios::binary ) ;
		if( !file.is_open() ) {
			throw FileNotOpenedError() ;
		}
		gen_file_ptr->push( file ); 
		return std::auto_ptr< std::ostream >( gen_file_ptr ) ;
  	}
}
