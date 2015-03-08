
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_FILE_UTILS_HPP
#define GENFILE_FILE_UTILS_HPP

#include <string>
#include <memory>
#include "genfile/Chromosome.hpp"

namespace genfile {
	struct CompressionType {
		CompressionType( char const* type ) ;
		CompressionType( std::string const& type ) ;
		CompressionType( CompressionType const& other ) ;
		CompressionType& operator=( CompressionType const& other ) ;
		bool operator==( CompressionType const& other ) const ;
	private:
		std::string m_type ;
		void check_type() const ;
	} ;

	CompressionType get_compression_type_indicated_by_filename( std::string const& filename ) ;
	Chromosome get_chromosome_indicated_by_filename( std::string const& filename ) ;

	std::pair< std::string, std::string > uniformise( std::string filename ) ;

	std::string get_gen_file_extension_if_present( std::string const& filename ) ;
	std::string strip_gen_file_extension_if_present( std::string const& filename ) ;
	std::pair< std::string, std::string > split_extension( std::string const& filename ) ;

	std::string replace_or_add_extension( std::string const& filename, std::string const& new_extension ) ;
	std::string create_temporary_filename() ;

	std::auto_ptr< std::istream > open_text_file_for_input( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::istream > open_text_file_for_input( std::string filename ) ;
	std::auto_ptr< std::ostream > open_text_file_for_output( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::ostream > open_text_file_for_output( std::string filename ) ;
	std::auto_ptr< std::istream > open_binary_file_for_input( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::istream > open_binary_file_for_input( std::string filename ) ;
	std::auto_ptr< std::ostream > open_binary_file_for_output( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::ostream > open_binary_file_for_output( std::string filename ) ;

	std::size_t count_lines_left_in_stream( std::istream& aStream, bool allow_empty_lines = false ) ;
}

#endif
