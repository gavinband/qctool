
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef DB_LOAD_KEY_VALUE_PAIRS_HPP
#define DB_LOAD_KEY_VALUE_PAIRS_HPP

#include <map>
#include <string>
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "genfile/Error.hpp"
#include "genfile/string_utils.hpp"

namespace db {
	template< typename Key, typename Value >
	std::map< Key, Value > load_key_value_pairs(
		Connection* connection,
		std::string const& table_or_view,
		std::string const& key_column,
		std::string const& value_column		
	)
	// Load key/value pairs from the given columns of the given table, returning the result as a std::map.
	// If a duplicate key is found, a KeyNotUniqueError is thrown.
	{
		typename std::map< Key, Value > result ;
		db::Connection::StatementPtr statement = connection->get_statement(
			"SELECT " + key_column + ", " + value_column + " FROM " + table_or_view + ";"
		) ;
	
		while( statement->step() ) {
			Key key = statement->get_column< Key >( 0 ) ;
			if( result.find( key ) != result.end() ) {
				throw genfile::DuplicateKeyError(
					"db::load_key_value_pairs()",
					connection->get_spec() + ":" + table_or_view + ", key=\"" + genfile::string_utils::to_string( key ) + "\""
				) ;
			}
			result.insert(
				std::make_pair(
			 		key,
					statement->get_column< Value >( 1 )
				)
			) ;
		}

		return result ;
	}
	
	template< typename Key, typename Value >
	std::map< Key, Value > load_key_value_pairs(
		Connection::UniquePtr& connection,
		std::string const& table_or_view,
		std::string const& key_column,
		std::string const& value_column		
	) {
		return load_key_value_pairs< Key, Value >(
			connection.get(),
			table_or_view,
			key_column,
			value_column
		) ;
	}
}

#endif
