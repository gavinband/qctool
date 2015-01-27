
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <sstream>
#include "test_case.hpp"
#include "sqlite3/sqlite3.h"
#include "db/Connection.hpp"
#include "db/SQLite3Connection.hpp"
#include "db/SQLite3Statement.hpp"
#include "db/Connection.hpp"
#include "db/load_key_value_pairs.hpp"

AUTO_TEST_CASE( test_load_key_value_pairs ) {
	std::cerr << "test_load_key_value_pairs()\n" ;
	db::Connection::UniquePtr connection( new db::SQLite3Connection( ":memory:" )) ;
	typedef db::Connection::StatementPtr StatementPtr ;
	
	// Construct a table
	{
		StatementPtr statement = connection->get_statement(
			"CREATE TABLE KeyValue ( key TEXT, value TEXT ) ;"
		) ;
		TEST_ASSERT( statement->step() == db::SQLite3Statement::Error::OK ) ;
	}
	
	{
		StatementPtr statement = connection->get_statement(
			"INSERT INTO KeyValue VALUES( ?1, ?2 ) ;"
		) ;
		
		for( std::size_t i = 0; i < 10; ++i ) {
			statement->reset() ;
			std::string key( i, 'x' ) ;
			std::string value( i*2, '*' ) ;
			statement->bind( 1, key ) ;
			statement->bind( 2, value ) ;
			TEST_ASSERT( statement->step() == db::SQLite3Statement::Error::OK ) ;
		}
	}
	typedef std::map< std::string, std::string > Map ;
	Map values = db::load_key_value_pairs< std::string, std::string >( connection, "KeyValue", "key", "value" ) ;
	TEST_ASSERT( values.size() == 10 ) ;
	for( std::size_t i = 0; i < 10; ++i ) {
		Map::const_iterator where = values.find( std::string( i, 'x' ) ) ;
		TEST_ASSERT( where != values.end() ) ;
		TEST_ASSERT( where->second == std::string( i*2, '*' )) ;
	}

	std::cerr << "success.\n" ;
}

AUTO_TEST_CASE( test_load_key_value_pairs_duplicates ) {
	std::cerr << "test_load_key_value_pairs_duplicates()\n" ;
	db::Connection::UniquePtr connection( new db::SQLite3Connection( ":memory:" )) ;
	typedef db::Connection::StatementPtr StatementPtr ;
	// Construct a table
	{
		StatementPtr statement = connection->get_statement(
			"CREATE TABLE KeyValue ( key TEXT, value TEXT ) ;"
		) ;
		TEST_ASSERT( statement->step() == db::SQLite3Statement::Error::OK ) ;
	}
	
	{
		StatementPtr statement = connection->get_statement(
			"INSERT INTO KeyValue VALUES( ?1, ?2 ) ;"
		) ;
		
		for( std::size_t i = 0; i < 10; ++i ) {
			statement->reset() ;
			statement->bind( 1, std::string( i, 'x' ) ) ;
			statement->bind( 2, std::string( i*2, '*' )) ;
			TEST_ASSERT( statement->step() == db::SQLite3Statement::Error::OK ) ;
		}
		
		statement->reset() ;
		statement->bind( 1, std::string( 1, 'x' )) ;
		TEST_ASSERT( statement->step() == db::SQLite3Statement::Error::OK ) ;
	}

	try {
		std::map< std::string, std::string > values = db::load_key_value_pairs< std::string, std::string >( connection, "KeyValue", "key", "value" ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::DuplicateKeyError const& e ) {
		// as expected.
	}
	std::cerr << "success.\n" ;
}


AUTO_TEST_MAIN {
	test_load_key_value_pairs() ;
	test_load_key_value_pairs_duplicates() ;
}

