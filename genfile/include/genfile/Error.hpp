
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_ERROR_HPP
#define GENFILE_ERROR_HPP

#include <exception>
#include <string>
#include <sstream>

namespace genfile {
	struct InputError: public std::exception 
	{
		virtual ~InputError() throw() {}

		InputError( std::string const& source, std::string const& message = "" ):
			m_source( source ),
			m_message( message )
		{}
		
		InputError( InputError const& other ):
			m_source( other.m_source ),
			m_message( other.m_message )
		{}

		const char* what() const throw() { return "genfile:InputError" ; }
		virtual std::string format_message() const ;
		std::string const& source() const { return m_source ; }
		std::string const& message() const { return m_message ; }
	private:
		std::string const m_source ;
		std::string const m_message ;
	} ;

	struct BadArgumentError: public InputError
	{
		BadArgumentError( std::string const& function, std::string const& arguments, std::string const& message = "" ):
			InputError( "arguments to " + function, message ),
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
		
		std::string format_message() const ;
		
		private:
			std::string const m_function ;
			std::string const m_arguments ;
	} ;

	struct MalformedInputError: public InputError
	{
		MalformedInputError( std::string const& source, std::string const& message, int line, int column = -1 ):
			InputError( source, message ),
			m_line( line ),
			m_column( column )
		{}
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

		std::string format_message() const ;
		
		int line() const { return m_line ; }
		bool has_column() const { return m_column >= 0 ; }
		int column() const { return m_column ; }

	private:
		int const m_line ;
		int const m_column ;
	} ;
	
	struct DuplicateEntryError: public InputError
	{
		~DuplicateEntryError() throw() {}

		DuplicateEntryError(
			 std::string const& caller, std::string const& source, std::string const& key, std::string const& variable, std::string const& message
		):
			InputError( source, message ),
			m_caller( caller ),
			m_key( key ),
			m_variable( variable )
		{}
		
		std::string const& key() const { return m_key ; }
		std::string const& variable() const { return m_variable ; }
		std::string const& caller() const { return m_caller ; }
		
		char const* what() const throw() { return "genfile::DuplicateEntryError" ; }
		std::string format_message() const ;
		
	private:
		std::string const& m_caller ;
		std::string const m_key ;
		std::string const m_variable ;
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

	struct FormatUnsupportedError: public InputError
	{
		FormatUnsupportedError( std::string const& source, std::string const& format ):
			InputError( source ),
			m_format( format )
		{}
		
		~FormatUnsupportedError() throw() {}
		
		char const* what() const throw() { return "genfile::FormatUnsupportedError" ; }

		std::string format_message() const ;
		
		std::string const& format() const {
			return m_format ;
		}
		
	private:
		std::string const m_format ;
	} ;

	struct DuplicateIndividualError: public MalformedInputError
	{
		DuplicateIndividualError( std::string const& source, std::string const& id_1, std::string const& message, int line, int column = -1 ):
			MalformedInputError( source, message, line, column ),
			m_id_1( id_1 )
		{}

		DuplicateIndividualError( DuplicateIndividualError const& other ):
			MalformedInputError( other ),
			m_id_1( other.m_id_1 )
		{}
		
		~DuplicateIndividualError() throw() {}
		
		char const* what() const throw() { return "genfile::DuplicateIndividualError" ; }
		std::string const& id_1() const { return m_id_1 ; }
		std::string format_message() const ;
	private:
		std::string const m_id_1 ;
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
			return MalformedInputError::format_message() + "\nThe file has two consecutive newlines.  Note that popular editors add an extra newline to the end of the file (which you can't see)." ;
		}

		char const* what() const throw() { return "genfile::FileHasTwoTrailingNewlinesError" ; }
	} ;

	struct MissingTrailingNewlineError: public MalformedInputError
	{
		MissingTrailingNewlineError( std::string const& source, int line ):
			MalformedInputError( source, line )
		{}
		
		MissingTrailingNewlineError( FileHasTwoTrailingNewlinesError const& other ):
			MalformedInputError( other )
		{}

		std::string format_message() const {
			return MalformedInputError::format_message() + "\nThe file should end in a newline." ;
		}

		char const* what() const throw() { return "genfile::MissingTrailingNewlineError" ; }
	} ;

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
		MismatchError( std::string const& caller, std::string const& source, std::string const& key1, std::string const& key2 ):
			m_caller( caller ),
			m_source( source ),
			m_key1( key1 ),
			m_key2( key2 )
		{}

		MismatchError( MismatchError const& other ):
			m_caller( other.m_caller ),
			m_source( other.m_source ),
			m_key1( other.m_key1 ),
			m_key2( other.m_key2 )
		{}
		
		~MismatchError() throw() {}
		char const* what() const throw() { return "genfile::MismatchError" ; }
		std::string const& source() const { return m_source ; } 
		std::string const& caller() const { return m_caller ; } 
		std::string const& key1() const { return m_key1 ; }
		std::string const& key2() const { return m_key2 ; }
	private:
		std::string const m_caller, m_source, m_key1, m_key2 ;
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
	
	struct OperationUnsupportedError: public std::exception
	{
		OperationUnsupportedError( std::string const& caller, std::string const& operation, std::string const& object ):
			m_caller( caller ),
			m_operation( operation ),
			m_object( object )
		{}

		OperationUnsupportedError( OperationUnsupportedError const& other ):
			m_operation( other.m_operation ),
			m_object( other.m_object )
		{}
		
		char const* what() const throw() { return "genfile::OperationUnsupportedError" ; }
		
		virtual ~OperationUnsupportedError() throw() {}
		
		std::string const& caller() const { return m_caller ;}
		std::string const& operation() const { return m_operation ;}
		std::string const& object() const { return m_object ;}
	private:
		std::string const m_caller ;
		std::string const m_operation ;
		std::string const m_object ;
	} ;

	struct OperationFailedError: public std::exception
	{
		OperationFailedError( std::string const& caller, std::string const& object, std::string const& operation ):
			m_caller( caller ),
			m_object( object ),
			m_operation( operation )
		{}

		OperationFailedError( OperationFailedError const& other ):
			m_caller( other.m_caller ),
			m_object( other.m_object ),
			m_operation( other.m_operation )
		{}
		
		virtual ~OperationFailedError() throw() {}

		char const* what() const throw() { return "genfile::OperationFailedError" ; }
		
		std::string const& caller() const { return m_caller ;}
		std::string const& operation() const { return m_operation ;}
		std::string const& object() const { return m_object ;}
		
		std::string get_message() const {
			return m_caller + " (" + m_object + "): operation " + m_operation + " failed." ;
		}

	private:
		std::string const m_caller ;
		std::string const m_object ;
		std::string const m_operation ;
	} ;
	


}

#endif
