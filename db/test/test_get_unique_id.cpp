#include <iostream>
#include <sstream>
#include "test_case.hpp"
#include "sqlite3/sqlite3.h"
#include "genfile/Error.hpp"
#include "db/Connection.hpp"
#include "db/SQLite3Connection.hpp"
#include "db/SQLite3Statement.hpp"
#include "db/Connection.hpp"
#include "db/get_unique_id.hpp"

AUTO_TEST_CASE( test_get_unique_id ) {
	std::cerr << "test_get_unique_id()\n" ;
	db::Connection::UniquePtr connection( new db::SQLite3Connection( ":memory:" )) ;
	typedef db::Connection::StatementPtr StatementPtr ;
	
	// Construct a table
	{
		StatementPtr statement = connection->get_statement(
			"CREATE TABLE Fields ( id INT, value INT ) ;"
		) ;
		statement->step() ;
	}
	
	{
		StatementPtr statement = connection->get_statement(
			"INSERT INTO Fields VALUES( ?1, ?2 ) ;"
		) ;
		
		for( std::size_t i = 0; i < 10; ++i ) {
			statement->reset() ;
			statement->bind( 1, i ) ;
			statement->bind( 2, i*2 ) ;
			TEST_ASSERT( statement->step() == db::SQLite3Statement::Error::OK ) ;
		}

		// make a duplicate key with value 10
		statement->reset() ;
		statement->bind( 1, 99 ) ;
		statement->bind( 2, 10 ) ;
		TEST_ASSERT( statement->step() == db::SQLite3Statement::Error::OK ) ;
	}
	
	for( std::size_t value = 0; value < 20; ++value ) {
		try {
			std::ostringstream ostr ;
			ostr << "WHERE value == " << value ;
			int id = db::get_unique_id( connection.get(), "Fields", "id", ostr.str() ) ;
			std::cerr << "value " << value << " found id " << id << ".\n" ;
			TEST_ASSERT( value % 2 == 0 ) ;
			TEST_ASSERT( id == value/2 ) ;
			TEST_ASSERT( id != 5 ) ;
		}
		catch( genfile::KeyNotFoundError const& e ) {
			TEST_ASSERT( value % 2 == 1 ) ;
			std::cerr << "value " << value << " found no ids.\n" ;
		}
		catch( genfile::DuplicateKeyError const& e ) {
			TEST_ASSERT( value == 10 ) ;
			std::cerr << "value " << value << " found duplicate ids.\n" ;
		}
	}
	std::cerr << "success.\n" ;
}

AUTO_TEST_MAIN {
	test_get_unique_id() ;
}