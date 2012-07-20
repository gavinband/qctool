
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef DB_MYSQL_STATEMENT_HPP
#define DB_MYSQL_STATEMENT_HPP

#include <cassert>
#include <string>
#include <exception>
#include "config.hpp"
#if HAVE_MYSQL
	#include "mysql.h"
#endif
#include "db/MySQLConnection.hpp"
#include "db/SQLStatement.hpp"

namespace db {
	class MySQLStatement: public SQLStatement
	{
	public:
		MySQLStatement( SQLite3Connection* connection, std::string const& SQL ) ;
		~MySQLStatement() ;
		bool step() ;
		bool empty() const ;
		
		std::size_t get_number_of_columns() const ;
		std::string get_name_of_column( std::size_t i ) const ;
		bool is_null( int column_id ) const ;
		std::string get_sql() const ;

		MySQLStatement& bind( std::size_t i, int32_t value ) ;
		MySQLStatement& bind( std::size_t i, uint32_t value ) ;
		MySQLStatement& bind( std::size_t i, int64_t value ) ;
		MySQLStatement& bind( std::size_t i, uint64_t value ) ;
		MySQLStatement& bind( std::size_t i, double value ) ;
		MySQLStatement& bind( std::size_t i, std::string const& value ) ;
		MySQLStatement& bind( std::size_t i, char const* buffer, char const* const end ) ;
		MySQLStatement& bind_NULL( std::size_t i ) ;
		MySQLStatement& reset() ;

	protected:
		
		int64_t get_column_int64( int column_id ) const ;
		int get_column_int( int column_id ) const ;
		double get_column_double( int column_id ) const ;
		std::string get_column_string( int column_id ) const ;
		char get_column_char( int column_id ) const ;
		
		int get_column_count() const ;

	private:
		std::string const m_sql ;
		std::auto_ptr< MYSQL_STMT > m_statement ;
		MySQLConnection* m_connection ;	
		bool m_executed ;
		bool m_have_results ;
	} ;
}

#endif
