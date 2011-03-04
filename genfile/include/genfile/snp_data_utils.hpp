#ifndef SNP_DATA_UTILS_HPP
#define SNP_DATA_UTILS_HPP

#define GENFILE_USE_FAST_PARSE_METHODS 1

#include <iostream>
#include <cassert>
#include <string>
#include <memory>
#include <vector>

#include "genfile/Chromosome.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
	
	enum CompressionType { e_NoCompression = 0, e_GzipCompression = 1 } ;

	bool filename_indicates_gen_format( std::string const& filename ) ;
	bool filename_indicates_bgen_format( std::string const& filename ) ;
	bool filename_indicates_gen_or_bgen_format( std::string const& filename ) ;
	CompressionType get_compression_type_indicated_by_filename( std::string const& filename ) ;
	Chromosome get_chromosome_indicated_by_filename( std::string const& filename ) ;

	std::string get_gen_file_extension_if_present( std::string const& filename ) ;
	std::string strip_gen_file_extension_if_present( std::string const& filename ) ;

	std::auto_ptr< std::istream > open_text_file_for_input( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::ostream > open_text_file_for_output( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::istream > open_binary_file_for_input( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::ostream > open_binary_file_for_output( std::string filename, CompressionType compression_type ) ;

	std::string create_temporary_filename() ;

	struct SNPDataError: public std::exception {
		SNPDataError( std::string const& filespec = "(unknown)" ): m_filespec( filespec ) {}
		~SNPDataError() throw() {}
		char const* what() const throw() { return "genfile::SNPDataError" ; }
		std::string const& get_filespec() const { return m_filespec ; }
	private:
		std::string const m_filespec ;
	} ;

	struct OperationFailedError: public std::exception { char const* what() const throw() { return "genfile::OperationFailedError" ; } } ;
	struct FileStructureInvalidError: public SNPDataError { char const* what() const throw() { return "genfile::FileStructureInvalidError" ; } } ;
	struct FileNotFoundError: public SNPDataError	
	{
		FileNotFoundError( std::string const& filespec ): m_filespec( filespec ) {}
		~FileNotFoundError() throw() {}
		char const* what() const throw() { return "genfile::FileNotFoundError" ; }
		std::string const& filespec() const { return m_filespec ; }
		private:
			std::string const m_filespec ;
	} ;
}

#endif
