#ifndef GENFILE_ERROR_HPP
#define GENFILE_ERROR_HPP

#include <exception>
#include <string>

namespace genfile {
	struct InputError: public std::exception 
	{
		InputError( std::string const& source ):
			m_source( source )
		{}
		virtual ~InputError() throw() {}
		std::string const& source() const { return m_source ; }
		const char* what() const throw() { return "genfile:InputError" ; }
	private:
		std::string const m_source ;
	} ;

	struct ResourceNotOpenedError: public InputError
	{
		ResourceNotOpenedError( std::string const& source ):
			InputError( source )
		{}
		ResourceNotOpenedError( ResourceNotOpenedError const& other ):
			InputError( other )
		{}
		~ResourceNotOpenedError() throw() {}

		char const* what() const throw() { return "genfile::ResourceNotOpenedError" ; }
	} ;

	struct MalformedInputError: public InputError
	{
		MalformedInputError():
			InputError( "(unknown)" ),
			m_line( -1 ),
			m_column( -1 )
		{}

		MalformedInputError( std::string const& source, int line, int column = -1 ):
			InputError( source ),
			m_line( line ),
			m_column( column )
		{}

		~MalformedInputError() throw() {}

		char const* what() const throw() { return "genfile::MalformedInputError" ; }
		int line() const { return m_line ; }
		bool has_column() const { return m_column >= 0 ; }
		int column() const { return m_column ; }

	private:
		int const m_line ;
		int const m_column ;
	} ;

	struct UnexpectedMissingValueError: public MalformedInputError
	{
		UnexpectedMissingValueError( std::string const& source, int line, int column = -1 ):
			MalformedInputError( source, line, column )
		{}

		char const* what() const throw() { return "genfile::UnexpectedMissingValueError" ; }
	} ;
	
	struct DuplicateKeyError: public MalformedInputError
	{
		DuplicateKeyError( std::string const& source, int line, int column = -1 ):
			MalformedInputError( source, line, column )
		{}

		char const* what() const throw() { return "genfile::DuplicateKeyError" ; }
	} ;

	struct FileHasTwoTrailingNewlinesError: public MalformedInputError
	{
		FileHasTwoTrailingNewlinesError( std::string const& source, int line ):
			MalformedInputError( source, line )
		{}

		char const* what() const throw() { return "genfile::FileHasTwoTrailingNewlinesError" ; }
	} ;

}

#endif
