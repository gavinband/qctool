
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <stdint.h>
#include <vector>
#include <string>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "qcdb/StorageOptions.hpp"
#include "qcdb/MultiVariantStorage.hpp"
#include "qcdb/FlatTableMultiVariantDBOutputter.hpp"
#include "qcdb/FlatFileMultiVariantOutputter.hpp"

namespace qcdb {
	MultiVariantStorage::UniquePtr MultiVariantStorage::create(
		std::string const& filename,
		std::size_t const number_of_key_fields,
		std::string const& analysis_name,
		std::string const& analysis_chunk,
		Metadata const& metadata,
		std::string const& compare_by,
		boost::optional< genfile::db::Connection::RowId > analysis_id	
	) {
		MultiVariantStorage::UniquePtr result ;
		std::vector< std::string > const file_spec = Storage::parse_filespec( filename ) ;
		if( file_spec[0] == "sqlite" ) {
			FlatTableMultiVariantDBOutputter::UniquePtr table_storage = FlatTableMultiVariantDBOutputter::create(
				file_spec[1],
				number_of_key_fields,
				analysis_name,
				analysis_chunk,
				metadata,
				compare_by,
				analysis_id
			) ;
			// Set table name, if specified
			if( file_spec.size() > 2 ) {
				table_storage->set_table_name( file_spec[2] ) ;
			}
			if( file_spec.size() > 3 ) {
				assert( file_spec[3] == "without rowid" ) ;
				table_storage->set_without_rowid() ;
			}
			result.reset( table_storage.release() ) ;
		} else if( file_spec[0] == "flat" ){
			FlatFileMultiVariantOutputter::UniquePtr file_storage = FlatFileMultiVariantOutputter::create(
				file_spec[1],
				number_of_key_fields,
				analysis_name,
				metadata
			) ;
			result.reset( file_storage.release() ) ;
		} else {
			assert(0) ;
		}

		return result ;
	}
}

