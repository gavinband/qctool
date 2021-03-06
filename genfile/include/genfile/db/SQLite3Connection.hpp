
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SQLITE3_CONNECTOR_HPP
#define GENFILE_SQLITE3_CONNECTOR_HPP

#include <cassert>
#include <string>
#include <exception>
#include "sqlite3/sqlite3.h"
#include "genfile/db/Connection.hpp"
#include "genfile/db/Transaction.hpp"
#include "genfile/db/Error.hpp"

extern "C" {
	void sqlite3_trace_callback( void* udp, const char* sql ) ;
	int sqlite3_busy_callback( void*, int number_of_tries ) ;
}

namespace genfile {
	namespace db {
		class SQLite3Connection: public Connection
		{
		public:
			// Open a connection to an on-disk database.
			// If the filename is ":memory:", open instead a connection
			// to a new private in-memory DB.
			SQLite3Connection( std::string const& filename, bool overwrite = true, std::string const& mode = "rw" ) ;
			virtual ~SQLite3Connection() ;

		public:

			StatementPtr get_statement( std::string const& SQL ) ;
			RowId get_last_insert_row_id() const ;

		public:
			sqlite3_stmt* prepare_sql( std::string const& SQL ) const ;
			int finalise_statement( sqlite3_stmt* statement ) ;
			int step_statement( sqlite3_stmt* statement ) ;
	
			std::string get_spec() const {
				return m_filename ;
			}

		private:
			struct Transaction: public db::Transaction {
				Transaction( SQLite3Connection& connection ) ;
				~Transaction() ;
			private:
				SQLite3Connection& m_connection ;
			} ;
		public:
			ScopedTransactionPtr open_transaction( double max_seconds_to_wait = 0.1 ) ;
		private:
			std::string m_filename ;
			sqlite3* m_db_connection ;
			bool m_managed ;

		private:
			void open_db_connection( std::string const& filename, bool overwrite, std::string const& mode = "rw" ) ;
			void close_db_connection_if_necessary() ;
			void finalise_prepared_statements() ;
		} ;
	}
}

#endif
