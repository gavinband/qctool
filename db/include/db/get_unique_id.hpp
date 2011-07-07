#ifndef DB_GET_UNIQUE_ID_HPP
#define DB_GET_UNIQUE_ID_HPP

#include "db/Connection.hpp"

namespace db {
	// Given a table, an "id" column of ints and a where clause,
	// return the id for the row matching the where clause.
	// throw a error::KeyError if there is not a unique matching row.
	int get_unique_id(
		Connection* connection,
		std::string const& table,
		std::string const& id_column,
		std::string const& where_clause
	) ;
}

#endif

