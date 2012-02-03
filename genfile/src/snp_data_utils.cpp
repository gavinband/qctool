#include <iostream>
#include <string>
#include <cstdio>
#include <vector>
#include <sstream>
#include <map>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include "genfile/snp_data_utils.hpp"
#include "genfile/Error.hpp"
#include "genfile/get_set.hpp"

namespace genfile {
	CompressionType::CompressionType( char const* type ):
		m_type( type )
	{
		check_type() ;
	}

	CompressionType::CompressionType( std::string const& type ):
		m_type( type )
	{
		check_type() ;
	}

	void CompressionType::check_type() const {
		if( m_type != "no_compression" && m_type != "gzip_compression" ) {
			throw BadArgumentError( "CompressionType::check_type()", "type==\"" + m_type + "\"." ) ;
		}
	}
	
	CompressionType::CompressionType( CompressionType const& other ):
		m_type( other.m_type )
	{}
	
	CompressionType& CompressionType::operator=( CompressionType const& other ) {
		m_type = other.m_type ;
		return *this ;
	}

	bool CompressionType::operator==( CompressionType const& other ) const {
		return m_type == other.m_type ;
	}
	
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
			return "gzip_compression" ;
		}
		else {
			return "no_compression" ;
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

	std::string replace_or_add_extension( std::string const& filename, std::string const& new_extension ) {
		std::string result ;
		std::size_t pos = filename.rfind( "." ) ;

		if( pos == std::string::npos ) {
			result = filename ;
		}
		else {
			result = filename.substr( 0, pos ) ;
		}
		if( new_extension.size() == 0 || new_extension[0] != '.' ) {
			result += "." ;
		}
		result += new_extension ;
		return result ;
	}

	std::string get_gen_file_extension_if_present( std::string const& filename ) {
		std::string recognised_extensions[5] = {
			".gen",
			".gen.gz",
			".bgen",
			".vcf",
			".vcf.gz"
		} ;
		
		for( std::size_t i = 0; i < 6u; ++i ) {
			if( filename.size() > recognised_extensions[i].size() && filename.substr( filename.size() - recognised_extensions[i].size(), recognised_extensions[i].size() ) == recognised_extensions[i] ) {
				return recognised_extensions[i] ;
			}
		}
		
		return "" ;
	}

	std::auto_ptr< std::istream > 
	open_text_file_for_input( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_istream > file_ptr( new boost::iostreams::filtering_istream ) ;
		if (compression_type == "gzip_compression") file_ptr->push(boost::iostreams::gzip_decompressor());
		boost::iostreams::file_source file( filename.c_str() ) ;
		if( !file.is_open() ) {
			throw ResourceNotOpenedError( filename ) ;
		}
		file_ptr->push( file ) ;
		return std::auto_ptr< std::istream >( file_ptr ) ;
  	}

	std::auto_ptr< std::istream > open_text_file_for_input( std::string filename ) {
		return open_text_file_for_input( filename, get_compression_type_indicated_by_filename( filename )) ;
	}

	std::auto_ptr< std::ostream > open_text_file_for_output( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_ostream > file_ptr( new boost::iostreams::filtering_ostream ) ;
		if( compression_type == "gzip_compression" ) {
			file_ptr->push( boost::iostreams::gzip_compressor() ) ;
		}

		if( filename == "-" ) {
			file_ptr->push( std::cout ) ;
		}
		else {
			boost::iostreams::file_sink file(filename.c_str()) ;
			if( !file.is_open() ) {
				throw ResourceNotOpenedError( filename ) ;
			}
			file_ptr->push( file ); 
		}
		return std::auto_ptr< std::ostream >( file_ptr ) ;
	}
	
	std::auto_ptr< std::ostream > open_text_file_for_output( std::string filename ) {
		return open_text_file_for_output( filename, get_compression_type_indicated_by_filename( filename )) ;
	}

	std::auto_ptr< std::istream > open_binary_file_for_input( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_istream > file_ptr( new boost::iostreams::filtering_istream ) ;
	    if (compression_type == "gzip_compression") file_ptr->push( boost::iostreams::gzip_decompressor() ) ;
		boost::iostreams::file_source file( filename.c_str(), std::ios::binary ) ;
		if( !file.is_open() ) {
			throw ResourceNotOpenedError( filename ) ;
		}
		file_ptr->push( file ) ;
		return std::auto_ptr< std::istream >( file_ptr ) ;
  	}

	std::auto_ptr< std::istream > open_binary_file_for_input( std::string filename ) {
		return open_binary_file_for_input( filename, get_compression_type_indicated_by_filename( filename )) ;
	}

	std::auto_ptr< std::ostream > open_binary_file_for_output( std::string filename, CompressionType compression_type ) {
		std::auto_ptr< boost::iostreams::filtering_ostream > file_ptr( new boost::iostreams::filtering_ostream ) ;
	    if (compression_type == "gzip_compression") file_ptr->push(boost::iostreams::gzip_compressor()) ;
		boost::iostreams::file_sink file(filename.c_str(), std::ios::binary ) ;
		if( !file.is_open() ) {
			throw ResourceNotOpenedError( filename ) ;
		}
		file_ptr->push( file ); 
		return std::auto_ptr< std::ostream >( file_ptr ) ;
  	}

	std::auto_ptr< std::ostream > open_binary_file_for_output( std::string filename ) {
		return open_binary_file_for_output( filename, get_compression_type_indicated_by_filename( filename )) ;
	}
	
	std::string create_temporary_filename() {
		return std::tmpnam( 0 ) ;
	}
	
	std::pair< std::string, std::string > uniformise( std::string filename ) {
		std::map< std::string, std::string > types  ;
		types[ ".bgen" ]							= "bgen" ;
		types[ ".gen" ] = types[ ".gen.gz" ] 		= "gen" ;
		types[ ".vcf" ] = types[ ".vcf.gz" ] 		= "vcf" ;
		types[ ".haps" ] = types[ ".haps.gz" ] 		= "impute_haplotypes" ;
		types[ ".phased" ] 							= "hapmap_haplotypes" ;

		for(
			std::map< std::string, std::string >::const_iterator i = types.begin();
			i != types.end();
			++i
		) {
			if(
				filename.size() >= i->first.size() &&
				filename.substr( filename.size()- i->first.size(), i->first.size() ) == i->first
			) {
				return std::make_pair( i->second, filename ) ;
			}
		}
		return std::make_pair( "", filename ) ;
	}
	
	std::size_t count_lines_left_in_stream( std::istream& aStream, bool allow_empty_lines ) {
		uint32_t number_of_lines = 0 ;
		std::vector< char > buffer( 1048576 ) ;
		do {
			aStream.read( &(buffer[0]), buffer.size() ) ;
			number_of_lines += std::count( buffer.begin(), buffer.begin() + aStream.gcount(), '\n' ) ;
			if( !allow_empty_lines ) {
				// A gen file can't contain a blank line.
				// Because popular editors (vim, nano, ..., but not emacs) typically add a trailing newline,
				// we might get in the situation where the GEN file has two trailing newlines thus messing
				// with our count.
				// Therefore we check here for the special case where what we've read ends in two newlines.
				if( (aStream.gcount() > 1) && (buffer[ aStream.gcount() - 1] == '\n') && (buffer[ aStream.gcount() - 2] == '\n') ) {
					throw FileHasTwoTrailingNewlinesError( "(stream)", number_of_lines ) ;
				}
			}
		}
		while( aStream ) ;

		// Most editors (vim, nano, but not emacs) automatically add a newline to the end of the file.
		// If the file has a trailing newline, we already have the correct count.
		// But if not, we've undercounted by one.
		if( aStream.gcount() > 0 ) {
			std::size_t pos = aStream.gcount() - 1 ;
			if( buffer[pos] != '\n' ) {
				++number_of_lines ;
			}
		}

		// We should have reached eof.
		// If so, report the data now.
		if( !aStream.eof() ) {
			throw MalformedInputError( "(stream)",  number_of_lines ) ;
		}
		return number_of_lines ;
	}
}

