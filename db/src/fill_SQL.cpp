#include <string>
#include "db/fill_SQL.hpp"

namespace db {
	std::string fill_SQL( std::string SQL ) {
		return SQL ;
	}

	std::string fill_SQL( std::string SQL, std::string const& value1 ) {
		std::size_t pos = SQL.find( "%s" ) ;
		assert( pos != std::string::npos ) ;
		SQL.replace( pos, 2, "\"" + value1 + "\"" ) ;
		return SQL ;
	}

	std::string fill_SQL( std::string SQL, char const* value1 ) {
		return fill_SQL( SQL, std::string( value1 )) ;
	}
}

