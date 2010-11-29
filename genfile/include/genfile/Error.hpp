#ifndef GENFILE_ERROR_HPP
#define GENFILE_ERROR_HPP

#include <exception>
#include <string>
#include <sstream>

namespace genfile {
	struct DuplicateKeyError: public std::exception
	{
		DuplicateKeyError( std::string const& source, std::string const& key ):
			m_source( source ),
			m_key( key )
		{}
		
		virtual ~DuplicateKeyError() throw() {}
	
		char const* what() const throw() { return "genfile::DuplicateKeyError" ; }
		std::string const& source() const { return m_source ; }
		std::string const& key() const { return m_key ; }
		
	private:
		std::string const m_source, m_key ;
	} ;

	struct DuplicateSNPError: public DuplicateKeyError
	{
		DuplicateSNPError( std::string const& source, std::string const& key ):
			DuplicateKeyError( source, key )
		{}
		char const* what() const throw() { return "genfile::DuplicateSNPError" ; }	
	} ;

	struct ColumnAlreadyExistsError: public DuplicateKeyError
	{
		ColumnAlreadyExistsError( std::string const& source, std::string const& key ):
			DuplicateKeyError( source, key )
		{}
		char const* what() const throw() { return "genfile::ColumnAlreadyExistsError" ; }	
	} ;
	
	struct MismatchError: public std::exception
		// Error thrown when two objects that should match (in some respect)
		// actually don't match.
	{
		MismatchError( std::string const& object1, std::string const& object2 ):
			m_object1( object1 ),
			m_object2( object2 )
		{}
		
		~MismatchError() throw() {}
		char const* what() const throw() { return "genfile::MismatchError" ; }
		std::string const& object1() const { return m_object1 ; }
		std::string const& object2() const { return m_object2 ; }
	private:
		
		std::string const m_object1, m_object2 ;
	} ;
	
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

	struct BadArgumentError: public InputError
	{
		BadArgumentError( std::string const& function, std::string const& arguments ):
			InputError( "arguments to " + function ),
			m_function( function ),
			m_arguments( arguments )
		{}
		
		~BadArgumentError() throw() {}
		
		char const* what() const throw() { return "genfile::BadArgumentError" ; }
		std::string const& function() const { return m_function ; }
		std::string const& arguments() const { return m_arguments ; }
		private:
			std::string const m_function ;
			std::string const m_arguments ;
	} ;

	struct KeyNotFoundError: public InputError
	{
		KeyNotFoundError( std::string const& key, std::string const& source ):
			InputError( source ),
			m_key( key)
		{}
		~KeyNotFoundError() throw() {}
		std::string const& key() const { return m_key ; }
		char const* what() const throw() { return "genfile::KeyNotFoundError" ; }
	private:
		
		std::string const m_key ;
	} ;
	
	struct ResourceNotOpenedError: public InputError
	{
		ResourceNotOpenedError( std::string const& source ):
			InputError( source )
		{}
		~ResourceNotOpenedError() throw() {}

		char const* what() const throw() { return "genfile::ResourceNotOpenedError" ; }
	} ;

	struct MalformedInputError: public InputError
	{
		MalformedInputError( std::string const& source, int line, int column = -1 ):
			InputError( source ),
			m_line( line ),
			m_column( column )
		{}

		~MalformedInputError() throw() {}

		char const* what() const throw() { return "genfile::MalformedInputError" ; }

		std::string format_message() const {
			std::ostringstream ostr ;
			ostr << "Source \"" << source() << "\" is malformed on line " << ( line() + 1 ) ;
			if( has_column() ) {
				ostr << ", column " + ( column() + 1 ) ;
			}
			ostr << "." ;
			return ostr.str() ;
		}
		
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
	
	
	struct DuplicateIndividualError: public MalformedInputError
	{
		DuplicateIndividualError( std::string const& source, std::string const& id_1, std::string const& id_2, int line, int column = -1 ):
			MalformedInputError( source, line, column ),
			m_id_1( id_1 ),
			m_id_2( id_2 )
		{}

		~DuplicateIndividualError() throw() {}
		
		char const* what() const throw() { return "genfile::DuplicateIndividualError" ; }
		std::string const& id_1() const { return m_id_1 ; }
		std::string const& id_2() const { return m_id_2 ; }
	private:
		std::string const m_id_1, m_id_2 ;
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
