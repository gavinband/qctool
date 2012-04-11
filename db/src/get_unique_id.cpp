
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "db/get_unique_id.hpp"
#include "genfile/Error.hpp"

namespace db {
	int get_unique_id(
		Connection* connection,
		std::string const& table,
		std::string const& id_column,
		std::string const& where_clause
	) {
		db::Connection::StatementPtr statement = connection->get_statement(
			"SELECT " + id_column + " FROM " + table + " " + where_clause + ";"
		) ;
		if( !statement->step() ) {
			throw genfile::DuplicateKeyError( "db::get_unique_id()", connection->get_spec() + ", <where clause> = \"" + where_clause + "\"" ) ;
		}
		int result = statement->get_column< int >( 0 ) ;
		if( statement->step() ) {
			throw genfile::DuplicateKeyError( "db::get_unique_id()", connection->get_spec() + ", <where clause> = \"" + where_clause + "\"" ) ;
		}
		return result ;
	}
}
