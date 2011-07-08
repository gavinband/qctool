#ifndef DB_SQL_STATEMENT_HPP
#define DB_SQL_STATEMENT_HPP

#include <cassert>
#include <string>
#include <exception>
#include "sqlite3/sqlite3.h"
#include "db/SQLite3Connection.hpp"

namespace db {
	
	struct SQLError: public std::exception
	{
		virtual ~SQLError() throw() {}
		char const* what() { return "db::SQLError" ; }
		virtual std::string description() const = 0 ;
	} ;
	
	class SQLStatement
	{
	public:
		typedef std::auto_ptr< SQLStatement > UniquePtr ;
	public:
		virtual ~SQLStatement() ;
		
		// Step to the next result row.  The entries can be accessed using get_column().
		// Return false if there is no next row, otherwise true.
		virtual bool step() = 0 ;
		// Return 
		virtual bool empty() const = 0 ;
		
		// Get the result for the given column.
		// Columns are 0-indexed.
		template< typename T > T get_column( int column_id ) const ;
		
		virtual std::size_t get_number_of_columns() const = 0 ;
		virtual std::string get_name_of_column( std::size_t i ) const = 0 ;

		// For parameterised queries, bind an integer value to the ith placeholder.
		// Placeholders are indexed starting from 1 on the left.
		// Named placeholders with the same name have the same index.
		virtual void bind( std::size_t i, int value ) const = 0 ;
		// For parameterised queries, bind a string value to the ith placeholder.
		// Placeholders are indexed starting from 1 on the left.
		// Named placeholders with the same name have the same index.
		// The string will be copied so the caller need not preserve it beyond the call site.
		virtual void bind( std::size_t i, std::string const& value ) const = 0 ;
		// For parameterised queries, bind a BLOB value (array of chars) to the ith placeholder.
		// Placeholders are indexed starting from 1 on the left.
		// Named placeholders with the same name have the same index.
		// The data will be copied and so the caller need not preserve it beyond the call site.
		virtual void bind( std::size_t i, char const* buffer, std::size_t n ) const = 0 ;

		// Reset the statement, ready to be re-executed.
		virtual void reset() const = 0 ;
		
		// Return the SQL this statement contains.
		virtual std::string get_sql() const = 0 ;
		
	protected:
		virtual int get_column_int( int column_id ) const = 0 ;
		virtual double get_column_double( int column_id ) const = 0 ;
		virtual std::string get_column_string( int column_id ) const = 0 ;
		virtual char get_column_char( int column_id ) const = 0 ;
	} ;
	
	template<> int SQLStatement::get_column< int >( int column_id ) const ;
	template<> double SQLStatement::get_column< double >( int column_id ) const ;
	template<> std::string SQLStatement::get_column< std::string >( int column_id ) const ;
}

#endif
