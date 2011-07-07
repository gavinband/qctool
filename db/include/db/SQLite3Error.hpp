#ifndef DB_SQLITE3_ERROR_HPP
#define DB_SQLITE3_ERROR_HPP

#include <string>
#include "db/SQLStatement.hpp"

namespace db {
	struct SQLite3Error: public std::exception
	{
		SQLite3Error( std::string const& caller, int error ) ;
		~SQLite3Error() throw() ;
		char const* what() const throw() { return "db::SQLite3Error" ; }

		int const& error() const { return m_error ; }
		std::string description() const ;
		
	private:
		std::string const m_caller ;
		int m_error ;
	} ;
}

#endif

