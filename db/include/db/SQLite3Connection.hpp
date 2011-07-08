#ifndef SQLITE3_CONNECTOR_HPP
#define SQLITE3_CONNECTOR_HPP

#include <cassert>
#include <string>
#include <exception>
#include "sqlite3/sqlite3.h"
#include "db/Connection.hpp"
#include "db/SQLite3Error.hpp"

namespace db {
	struct SQLite3OpenDBError: public SQLite3Error
	{
		SQLite3OpenDBError( std::string const& caller, int error_code ): SQLite3Error( caller, error_code ) {}
		char const* what() const throw() { return "db::SQLite3OpenDBError" ; }
	} ;

	struct SQLite3PrepareStatementError: public SQLite3Error
	{
		SQLite3PrepareStatementError( std::string const& caller, int error_code, std::string SQL ): SQLite3Error( caller, error_code ), m_SQL( SQL ) {}
		~SQLite3PrepareStatementError() throw() {}
		char const* what() const throw() { return "db::SQLite3PrepareStatementError" ; }
		std::string const& get_SQL() const { return m_SQL ; }
	private:
		std::string const m_SQL ;
	} ;

	struct SQLite3StepStatementError: public SQLite3Error
	{
		SQLite3StepStatementError( std::string const& caller, int error_code ): SQLite3Error( caller, error_code ) {}
		char const* what() const throw() { return "db::SQLite3StepStatementError" ; }
	} ;

	class SQLite3Connection: public Connection
	{
	public:
		// Open a connection to an on-disk database.
		// If the filename is ":memory:", open instead a connection
		// to a new private in-memory DB.
		SQLite3Connection( std::string const& filename ):
			m_filename( filename ),
			m_db_connection(0),
			m_managed( true )
		{
			open_db_connection( filename ) ;
		}

		virtual ~SQLite3Connection() {
			close_db_connection_if_necessary() ;
		}

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
		std::string m_filename ;
		sqlite3* m_db_connection ;
		bool m_managed ;

	private:
		void open_db_connection( std::string const& filename ) ;
		void close_db_connection_if_necessary() ;
		void finalise_prepared_statements() ;
	} ;
}

#endif
