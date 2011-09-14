#include <string>
#include "db/Connection.hpp"
#include "DataStore.hpp"
#include "VCDBDataStore.hpp"

DataStore::UniquePtr DataStore::create( std::string const& filename ) {
	DataStore::UniquePtr result ;
	db::Connection::UniquePtr connection( db::Connection::create( filename )) ;
	result.reset(
		new VCDBDataStore( connection )
	) ;
	return result ;
}
