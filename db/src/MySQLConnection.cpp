
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
#if 0

#include <iostream>
#include <stdint.h>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/thread/thread_time.hpp>
#include <boost/thread/thread.hpp>
#include "db/MySQLConnection.hpp"
// #include "db/MySQLStatement.hpp"
#include "db/Error.hpp"

namespace db {
	struct MySQLStatement ;

	MySQLConnection::MySQLConnection( std::string const& filename, bool overwrite ):
		m_filename( filename ),
		m_db_connection(0),
		m_managed( true )
	{
		open_db_connection( filename, overwrite ) ;
	}

	MySQLConnection::~MySQLConnection() {
		close_db_connection_if_necessary() ;
	}

	
	MySQLConnection::StatementPtr MySQLConnection::get_statement( std::string const& SQL ) {
		return StatementPtr( new MySQLStatement( this, SQL ) ) ;
	}	

	Connection::RowId MySQLConnection::get_last_insert_row_id() const {
		uint64_t result = sqlite3_last_insert_rowid( m_db_connection ) ;
		return result ;
	}

	void MySQLConnection::open_db_connection( std::string const& connection_string, bool overwrite ) {
		std::map< std::string, std::string > bits = parse_connection_string( connection_string ) ;
		open_db_connection( connection_spec, overwrite ) ;
	}


	void MySQLConnection::open_db_connection( std::map< std::string, std::string > const& connection_spec, bool overwrite ) {
		int flags = SQLITE_OPEN_READWRITE ;
		if( overwrite ) {
			flags |= SQLITE_OPEN_CREATE ;
		}
		int code = sqlite3_open_v2( filename.c_str(), &m_db_connection, flags, NULL ) ;
		if( code != SQLITE_OK ) {
			throw ConnectionError( "MySQLConnection::open_db_connection()", get_spec(), code ) ;
		}
		// We add a busy handler.  This makes the database more robust by retrying failed transactions.
		sqlite3_busy_handler( m_db_connection, &sqlite3_busy_callback, NULL ) ;
		// Uncomment the next line to trace SQL statements executed.
		//sqlite3_trace( m_db_connection, &sqlite3_trace_callback, NULL ) ;
	}

	void MySQLConnection::close_db_connection_if_necessary() {
		if( m_managed && m_db_connection != 0 ) {
			// According to SQLite docs, we must finalise any prepared statements
			// before we can close the db connection.
			finalise_prepared_statements() ;
			sqlite3_close( m_db_connection ) ;
			m_db_connection = 0 ;
		}
	}

	void MySQLConnection::finalise_prepared_statements() {
		assert( m_db_connection != 0 ) ;
#if SQLITE_VERSION_NUMBER > 3006000
		sqlite3_stmt* statement ;
		while(( statement = sqlite3_next_stmt( m_db_connection, 0)) != 0 ) {
			finalise_statement( statement ) ;
		}
#endif
	}
	
	MySQLConnection::ScopedTransactionPtr MySQLConnection::open_transaction( double max_seconds_to_wait ) {		
		ScopedTransactionPtr transaction ;
		for( std::size_t i = 0; i < max_seconds_to_wait * 100; ++i ) {
			try {
				transaction.reset( new MySQLConnection::Transaction( *this ) ) ;
				break ;
			}
			catch( db::StatementStepError const& e ) {
				boost::this_thread::sleep( boost::posix_time::milliseconds( 10 ) ) ;
			}
		}
		return transaction ;
	}
	
	MySQLConnection::Transaction::Transaction( MySQLConnection& connection ):
		m_connection( connection )
	{
		m_connection.run_statement( "BEGIN IMMEDIATE TRANSACTION" ) ;
	}
	
	MySQLConnection::Transaction::~Transaction() {
		m_connection.run_statement( "COMMIT" ) ;
	}
}

#endif
