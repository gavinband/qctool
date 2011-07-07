#include <cassert>
#include <string>
#include "sqlite3/sqlite3.h"
#include "db/SQLStatement.hpp"

namespace db {
	
	SQLStatement::~SQLStatement() {}
	
	template<>
	int SQLStatement::get_column< int >( int column_id ) const {
		return this->get_column_int( column_id ) ;
	}

	template<>
	double SQLStatement::get_column< double >( int column_id ) const {
		return this->get_column_double( column_id ) ;
	}

	template<>
	std::string SQLStatement::get_column< std::string >( int column_id ) const {
		return this->get_column_string( column_id ) ;
	}

	template<>
	char SQLStatement::get_column< char >( int column_id ) const {
		return this->get_column_char( column_id ) ;
	}
}
