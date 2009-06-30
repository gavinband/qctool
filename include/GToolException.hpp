#ifndef __GTOOL_EXCEPTION__
#define __GTOOL_EXCEPTION__

// base class of gtool exceptions.

#include <iostream>
#include <exception>

struct GToolException: public std::exception
{
	enum { BUFFER_SIZE = 2000 } ;
	
	public:
		GToolException( std::string const& msg = "" ) {
			set_message( msg ) ;
		}

	protected:
		void set_message( std::string const& msg ) const throw() {
			strncpy( m_buffer, msg.c_str(), std::min( msg.size(), static_cast< std::size_t > ( BUFFER_SIZE ))) ;
			m_buffer[ BUFFER_SIZE - 1 ] = '\0' ;
			}

	public:
		// From std::exception.
		char const* what() const throw() { return m_buffer ; }


	private:
		static char m_buffer[ BUFFER_SIZE ] ;
	
		friend std::ostream& operator<<( std::ostream&, GToolException const& ) ;
} ;


struct BadFileFormatException: public GToolException
{
	BadFileFormatException( std::string const& msg )
	: GToolException( msg )
	{}
} ;

struct BadRowFormatException: public GToolException
{
	BadRowFormatException( std::string const& msg )
	: GToolException( msg )
	{}
} ;

#endif

