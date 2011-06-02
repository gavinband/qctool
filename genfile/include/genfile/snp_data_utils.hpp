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
#include "genfile/FileUtils.hpp"

namespace genfile {
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
