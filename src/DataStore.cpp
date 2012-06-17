
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include "genfile/Error.hpp"
#include "db/Connection.hpp"
#include "DataStore.hpp"
#include "VCDBDataStore.hpp"

DataStore::UniquePtr DataStore::create( std::string const& spec ) {
	DataStore::UniquePtr result ;
	if( spec.compare( 0, 10, "sqlite3://" ) == 0 ) {
		db::Connection::UniquePtr connection( db::Connection::create( spec.substr( 10, spec.size() ))) ;
		result.reset(
			new VCDBDataStore( connection )
		) ;
	} else {
		throw genfile::BadArgumentError( "DataStore::create()", "spec=\"" + spec + "\"" ) ;
	}
	return result ;
}
