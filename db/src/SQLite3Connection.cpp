#include <stdint.h>
#include "db/SQLite3Connection.hpp"
#include "db/SQLite3Statement.hpp"

namespace db {
	SQLite3Connection::StatementPtr SQLite3Connection::get_statement( std::string const& SQL ) {
		return StatementPtr( new SQLite3Statement( this, SQL ) ) ;
	}	

	Connection::RowId SQLite3Connection::get_last_insert_row_id() const {
		uint64_t result = sqlite3_last_insert_rowid( m_db_connection ) ;
		return result ;
	}

	sqlite3_stmt* SQLite3Connection::prepare_sql( std::string const& SQL ) const {
		assert( m_db_connection != 0 ) ;
		sqlite3_stmt* statement ;
		int code = sqlite3_prepare_v2(
			m_db_connection,
			SQL.c_str(),
			static_cast< int >( SQL.size()+1 ), // +1 accounts for null terminating byte.
			&statement,
			0 // ignore pzTail
		) ;
		if( code != SQLITE_OK ) {
			throw SQLite3PrepareStatementError( "SQLite3Connection::prepare_sql()", code, SQL ) ;
		}
		// SQLite might return 0 if the SQL consisted of comments and whitespace only
		// but I want to treat this as a programmer error.
		assert( statement != 0 ) ;
		return statement ;
	}

	int SQLite3Connection::finalise_statement( sqlite3_stmt* statement ) {
		assert( statement != 0 ) ;
		return sqlite3_finalize( statement ) ;
	}

	int SQLite3Connection::step_statement( sqlite3_stmt* statement ) {
		int code = sqlite3_step( statement ) ;
		if( code != SQLITE_ROW && code != SQLITE_DONE ) {
			throw SQLite3StepStatementError( "SQLite3Connection::step_statement()", code ) ;
		}
		return (code == SQLITE_ROW);
	}
	
	void SQLite3Connection::open_db_connection( std::string const& filename ) {
		int code = sqlite3_open( filename.c_str(), &m_db_connection ) ;
		if( code != SQLITE_OK ) {
			throw SQLite3OpenDBError( "SQLite3Connection::open_db_connection()", code ) ;
		}
	}

	void SQLite3Connection::close_db_connection_if_necessary() {
		if( m_managed && m_db_connection != 0 ) {
			// According to SQLite docs, we must finalise any prepared statements
			// before we can close the db connection.
			finalise_prepared_statements() ;
			sqlite3_close( m_db_connection ) ;
			m_db_connection = 0 ;
		}
	}

	void SQLite3Connection::finalise_prepared_statements() {
		assert( m_db_connection != 0 ) ;
#if SQLITE_VERSION_NUMBER > 3006000
		sqlite3_stmt* statement ;
		while(( statement = sqlite3_next_stmt( m_db_connection, 0)) != 0 ) {
			finalise_statement( statement ) ;
		}
#endif
	}
}
