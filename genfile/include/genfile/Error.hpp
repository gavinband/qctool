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

		DuplicateKeyError( DuplicateKeyError const& other ):
			m_source( other.m_source ),
			m_key( other.m_key )
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
		DuplicateSNPError( DuplicateSNPError const& other ):
			DuplicateKeyError( other )
		{}

		char const* what() const throw() { return "genfile::DuplicateSNPError" ; }	
	} ;

	struct ColumnAlreadyExistsError: public DuplicateKeyError
	{
		ColumnAlreadyExistsError( std::string const& source, std::string const& key ):
			DuplicateKeyError( source, key )
		{}
		
		ColumnAlreadyExistsError( ColumnAlreadyExistsError const& other ):
			DuplicateKeyError( other )
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

		MismatchError( MismatchError const& other ):
			m_object1( other.m_object1 ),
			m_object2( other.m_object2 )
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
		
		InputError( InputError const& other ):
			m_source( other.m_source )
		{}

		virtual std::string format_message() const { return "An unknown input error occurred in source \"" + source() + "\"." ; }

		virtual ~InputError() throw() {}
		std::string const& source() const { return m_source ; }
		const char* what() const throw() { return "genfile:InputError" ; }
	private:
		std::string const m_source ;
	} ;

	struct OutputError: public std::exception 
	{
		OutputError( std::string const& source ):
			m_source( source )
		{}
		
		OutputError( OutputError const& other ):
			m_source( other.m_source )
		{}

		virtual std::string format_message() const { return "An unknown output error occurred in source \"" + source() + "\"." ; }

		virtual ~OutputError() throw() {}
		std::string const& source() const { return m_source ; }
		const char* what() const throw() { return "genfile:OutputError" ; }
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

		BadArgumentError( BadArgumentError const& other ):
			InputError( other ),
			m_function( other.m_function ),
			m_arguments( other.m_arguments )
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
		KeyNotFoundError( std::string const& source, std::string const& key ):
			InputError( source ),
			m_key( key)
		{}

		KeyNotFoundError( KeyNotFoundError const& other ):
			InputError( other ),
			m_key( other.m_key )
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
		ResourceNotOpenedError( ResourceNotOpenedError const& other ):
			InputError( other )
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

		MalformedInputError( MalformedInputError const& other ):
			InputError( other ),
			m_line( other.m_line ),
			m_column( other.m_column )
		{}

		~MalformedInputError() throw() {}

		char const* what() const throw() { return "genfile::MalformedInputError" ; }

		virtual std::string format_message() const {
			std::ostringstream ostr ;
			ostr << "Source \"" << source() << "\" is malformed on line " << ( line() + 1 ) ;
			if( has_column() ) {
				ostr << ", column " << ( column() + 1 ) ;
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

	struct FormatUnsupportedError: public InputError
	{
		FormatUnsupportedError( std::string const& source, std::string const& format ):
			InputError( source ),
			m_format( format )
		{}
		
		~FormatUnsupportedError() throw() {}
		
		char const* what() const throw() { return "genfile::FormatUnsupportedError" ; }

		virtual std::string format_message() const {
			return "Source \"" + source() + "\" is in format \"" + format() + "\", which is unsupported." ;
		}
		
		std::string const& format() const {
			return m_format ;
		}
		
	private:
		std::string const m_format ;
	} ;

	struct UnexpectedMissingValueError: public MalformedInputError
	{
		UnexpectedMissingValueError( std::string const& source, int line, int column = -1 ):
			MalformedInputError( source, line, column )
		{}
		
		UnexpectedMissingValueError( UnexpectedMissingValueError const& other ):
			MalformedInputError( other )
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

		DuplicateIndividualError( DuplicateIndividualError const& other ):
			MalformedInputError( other ),
			m_id_1( other.m_id_1 ),
			m_id_2( other.m_id_2 )
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
		
		FileHasTwoTrailingNewlinesError( FileHasTwoTrailingNewlinesError const& other ):
			MalformedInputError( other )
		{}

		std::string format_message() const {
			return MalformedInputError::format_message() + "\nThe file has two consecutive newlines.  Note that popular editors add an extra newline to the end of the file (which you can't see).\n" ;
		}

		char const* what() const throw() { return "genfile::FileHasTwoTrailingNewlinesError" ; }
	} ;
}

#endif
