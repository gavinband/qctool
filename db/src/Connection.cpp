#include "db/Connection.hpp"
#include "db/SQLite3Connection.hpp"

namespace db {
	Connection::UniquePtr Connection::create( std::string const& filename ) {
		return Connection::UniquePtr( new SQLite3Connection( filename )) ;
	}

	void Connection::run_statement( std::string const& SQL ) {
		get_statement( SQL )->step() ;
	}
}
