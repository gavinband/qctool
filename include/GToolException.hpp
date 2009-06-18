#ifndef __GTOOL_EXCEPTION__
#define __GTOOL_EXCEPTION__

// base class of gtool exceptions.

#include <iostream>
#include <exception>

struct GToolException: public std::exception
{
	public:
		GToolException( std::string const& msg ) {
			strncpy( m_buffer, msg.c_str(), 2000 ) ;
		}

		char const * const message() const { return m_buffer; }

	private:
		static char m_buffer[2000] ;
	
		friend std::ostream& operator<<( std::ostream&, GToolException const& ) ;
};

struct BadRowFormatException: public GToolException
{
	BadRowFormatException( std::string const& msg )
	: GToolException( msg )
	{}
} ;

#endif

