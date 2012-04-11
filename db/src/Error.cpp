
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include "db/SQLite3Statement.hpp"
#include "db/Error.hpp"

namespace db {
	Error::Error( std::string const& caller, std::string const& db_spec, int error ):
		m_spec( db_spec ),
		m_error( error )
	{
		assert( m_error != SQLite3Statement::Error::OK ) ;
	}
	
	Error::~Error() throw() {}
	
	std::string Error::description() const {
		typedef db::SQLite3Statement::Error Error ;
		switch( m_error ) {
			case Error::ERROR: 		return "SQL error or missing database"; break ;
			case Error::INTERNAL: 	return "Internal logic error in SQLite"; break ;
			case Error::PERM: 		return "Access permission denied"; break ;
			case Error::ABORT: 		return "Callback routine requested an abort"; break ;
			case Error::BUSY: 		return "The database file is locked"; break ;
			case Error::LOCKED: 	return "A table in the database is locked"; break ;
			case Error::NOMEM: 		return "A malloc() failed"; break ;
			case Error::READONLY: 	return "Attempt to write a readonly database"; break ;
			case Error::INTERRUPT: 	return "Operation terminated by sqlite3_interrupt()"; break ;
			case Error::IOERR: 		return "Some kind of disk I/O error occurred"; break ;
			case Error::CORRUPT: 	return "The database disk image is malformed"; break ;
			case Error::NOTFOUND: 	return "NOT USED. Table or record not found"; break ;
			case Error::FULL: 		return "Insertion failed because database is full"; break ;
			case Error::CANTOPEN: 	return "Unable to open the database file"; break ;
			case Error::PROTOCOL: 	return "Database lock protocol error"; break ;
			case Error::EMPTY: 		return "Database is empty"; break ;
			case Error::SCHEMA: 	return "The database schema changed"; break ;
			case Error::TOOBIG: 	return "String or BLOB exceeds size limit"; break ;
			case Error::CONSTRAINT: return "Abort due to constraint violation"; break ;
			case Error::MISMATCH: 	return "Data type mismatch"; break ;
			case Error::MISUSE: 	return "Library used incorrectly"; break ;
			case Error::NOLFS: 		return "Uses OS features not supported on host"; break ;
			case Error::AUTH: 		return "Authorization denied"; break ;
			case Error::FORMAT: 	return "Auxiliary database format error"; break ;
			case Error::RANGE: 		return "2nd parameter to sqlite3_bind out of range"; break ;
			case Error::NOTADB: 	return "File opened that is not a database file"; break ;
			case Error::ROW: 		return "sqlite3_step() has another row ready"; break ;
			case Error::DONE: 		return "sqlite3_step() has finished executing"; break ;
			default: 				assert(0) ;
		}
		return "(unknown)" ;
	}
}

