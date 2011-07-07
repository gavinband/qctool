#ifndef DB_CONNECTION_HPP
#define DB_CONNECTION_HPP

#include <memory>
#include <string>

namespace db {

	class SQLStatement ;
	
	class Connection
		// Base class for classes representing a connection to a database.
		// The only supported operation is getting a query representing some SQL.
	{
	public:
		typedef std::auto_ptr< Connection > UniquePtr ;
		typedef std::auto_ptr< SQLStatement > StatementPtr ;
		
		static UniquePtr create( std::string const& filename ) ;
		
		virtual ~Connection() {}

		virtual StatementPtr get_statement( std::string const& SQL ) = 0 ;
		virtual std::string get_spec() const = 0 ;
	} ;	
}

#endif
