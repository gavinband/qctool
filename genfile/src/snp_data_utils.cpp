#include <iostream>
#include <string>
#include <cstdio>
#include <vector>
#include <sstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include "genfile/snp_data_utils.hpp"
#include "genfile/Error.hpp"
#include "genfile/get_set.hpp"

namespace genfile {
	bool filename_indicates_gen_format( std::string const& filename ) {
		return ( filename.find( ".gen") != std::string::npos ) ;
	}
	bool filename_indicates_bgen_format( std::string const& filename ) {
		return ( filename.find( ".bgen") != std::string::npos ) ;
	}
	bool filename_indicates_gen_or_bgen_format( std::string const& filename ) {
		return filename_indicates_gen_format( filename ) || filename_indicates_bgen_format( filename ) ;
	}

	CompressionType get_compression_type_indicated_by_filename( std::string const& filename ) {
		if( filename.find( ".gz") != std::string::npos ) {
			return e_GzipCompression ;
		}
		else {
			return e_NoCompression ;
		}
	}

	Chromosome get_chromosome_indicated_by_filename( std::string const& filename ) {
		std::vector< Chromosome > indicated_chromosomes ;
		std::ostringstream ostr ;
		for( unsigned char i = 0; i < 23; ++i ) {
			ostr.str( "" ) ;
			ostr << Chromosome( i ) ;
			if( filename.find( ostr.str() ) != std::string::npos ) {
				indicated_chromosomes.push_back( Chromosome( i )) ;
			}
		}
		ostr.str( "" ) ;
		ostr << XYPseudoAutosomalDNA ;
		if( filename.find( ostr.str() ) != std::string::npos ) {
			indicated_chromosomes.push_back( XYPseudoAutosomalDNA ) ;
		}
		ostr.str( "" ) ;
		ostr << MitochondrialDNA ;
		if( filename.find( ostr.str() ) != std::string::npos ) {
			indicated_chromosomes.push_back( MitochondrialDNA ) ;
		}

		if( indicated_chromosomes.size() != 1 ) {
			return UnidentifiedChromosome ;
		}

		return indicated_chromosomes[0] ;
	}

	std::string strip_gen_file_extension_if_present( std::string const& filename ) {
		std::string extension = get_gen_file_extension_if_present( filename ) ;
		return filename.substr( 0, filename.size() - extension.size() ) ;
	}

	std::string get_gen_file_extension_if_present( std::string const& filename ) {
		std::string recognised_extensions[4] = {
			".gen",
			".gen.gz",
			".bgen",
			".bgen.gz"
		} ;
		
		for( std::size_t i = 0; i < 4u; ++i ) {
			if( filename.substr( filename.size() - recognised_extensions[i].size(), recognised_extensions[i].size() ) == recognised_extensions[i] ) {
				return recognised_extensions[i] ;
			}
		}

		return "" ;
	}

	std::auto_ptr< std::istream > open_text_file_for_input( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_istream > gen_file_ptr( new boost::iostreams::filtering_istream ) ;
	    if (compression_type == e_GzipCompression) gen_file_ptr->push(boost::iostreams::gzip_decompressor());
		boost::iostreams::file_source file( filename.c_str() ) ;
		if( !file.is_open() ) {
			throw ResourceNotOpenedError( filename ) ;
		}
		gen_file_ptr->push( file ) ;
		return std::auto_ptr< std::istream >( gen_file_ptr ) ;
  	}

	std::auto_ptr< std::ostream > open_text_file_for_output( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_ostream > gen_file_ptr( new boost::iostreams::filtering_ostream ) ;
	    if (compression_type == e_GzipCompression) gen_file_ptr->push(boost::iostreams::gzip_compressor());
		boost::iostreams::file_sink file(filename.c_str()) ;
		if( !file.is_open() ) {
			throw ResourceNotOpenedError( filename ) ;
		}
		gen_file_ptr->push( file ); 
		return std::auto_ptr< std::ostream >( gen_file_ptr ) ;
  	}

	std::auto_ptr< std::istream > open_binary_file_for_input( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_istream > gen_file_ptr( new boost::iostreams::filtering_istream ) ;
	    if (compression_type == e_GzipCompression) gen_file_ptr->push( boost::iostreams::gzip_decompressor() ) ;
		boost::iostreams::file_source file( filename.c_str(), std::ios::binary ) ;
		if( !file.is_open() ) {
			throw ResourceNotOpenedError( filename ) ;
		}
		gen_file_ptr->push( file ) ;
		return std::auto_ptr< std::istream >( gen_file_ptr ) ;
  	}

	std::auto_ptr< std::ostream > open_binary_file_for_output( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_ostream > gen_file_ptr( new boost::iostreams::filtering_ostream ) ;
	    if (compression_type == e_GzipCompression) gen_file_ptr->push(boost::iostreams::gzip_compressor()) ;
		boost::iostreams::file_sink file(filename.c_str(), std::ios::binary ) ;
		if( !file.is_open() ) {
			throw ResourceNotOpenedError( filename ) ;
		}
		gen_file_ptr->push( file ); 
		return std::auto_ptr< std::ostream >( gen_file_ptr ) ;
  	}

	std::string create_temporary_filename() {
		return std::tmpnam( 0 ) ;
	}
	
	std::pair< std::string, std::string > uniformise( std::string filename ) {
		if( filename.find( ".gen") != std::string::npos ) {
			return std::make_pair( "gen", filename ) ;
		}
		else if( filename.find( ".bgen" ) != std::string::npos ) {
			return std::make_pair( "bgen", filename ) ;
		}
		else if( filename.find( ".vcf" ) != std::string::npos ) {
			return std::make_pair( "vcf", filename ) ;
		}
		else {
			return std::make_pair( "", filename ) ;
		}
	}
}
