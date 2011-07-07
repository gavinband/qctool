#include "db/Connection.hpp"
#include "db/SQLite3Connection.hpp"

namespace db {
	Connection::UniquePtr Connection::create( std::string const& filename ) {
		return Connection::UniquePtr( new SQLite3Connection( filename )) ;
	}
}
