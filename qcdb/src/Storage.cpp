
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
#include "qcdb/Storage.hpp"
#include "qcdb/FlatTableDBOutputter.hpp"
#include "qcdb/FlatFileOutputter.hpp"

namespace qcdb {
	std::vector< std::string > Storage::parse_filespec( std::string spec ) {
		std::vector< std::string > result(2) ;
		result[0] = "flat" ;
	
		// First get rid of leading sqlite specifier, if any.
		if( spec.size() >= 9 && spec.substr( 0, 9 ) == "sqlite://" ) {
			result[0] = "sqlite" ;
			spec = spec.substr( 9, spec.size() ) ;
		}

		std::vector< std::string > elts = genfile::string_utils::split( spec, ":" ) ;
		assert( elts.size() > 0 ) ;

		if( elts[0].size() >= 7 && elts[0].substr( elts[0].size() - 7, 7 ) == ".sqlite" ) {
			result[0] = "sqlite" ;
		}
		result[1] = elts[0] ;

		if( elts.size() > 2 ) {
			throw genfile::BadArgumentError(
				"parse_filespec()",
				"spec=\"" + spec + "\"",
				"Expected format for filespec is <filename> or sqlite://<filename>[:<tablename>]."
			) ;
		}
		if( elts.size() == 2 ) {
			if( result[0] != "sqlite" ) {
				throw genfile::BadArgumentError(
					"parse_filespec()",
					"spec=\"" + spec + "\"",
					"Table spec is not expected unless a sqlite file is specified (i.e. sqlite://<filename>:<tablename>)"
				) ;
			}
			std::vector< std::string > table_elts = genfile::string_utils::split( elts[1], "+" ) ;
			if( table_elts.size() == 2 ) {
				result.push_back( table_elts[0] ) ;
				result.push_back( table_elts[1] ) ;
			} else {
				result.push_back( elts[1] ) ;
			}
		}
	
		return( result ) ;
	}

	Storage::UniquePtr Storage::create(
		std::string const& filename,
		std::string const& analysis_name,
		std::string const& analysis_chunk,
		Metadata const& metadata,
		std::string const& compare_by,
		boost::optional< genfile::db::Connection::RowId > analysis_id
	) {
		Storage::UniquePtr result ;
		std::vector< std::string > const file_spec = Storage::parse_filespec( filename ) ;
		if( file_spec[0] == "sqlite" ) {
			FlatTableDBOutputter::UniquePtr table_storage = FlatTableDBOutputter::create(
				file_spec[1],
				analysis_name,
				analysis_chunk,
				metadata,
				compare_by,
				analysis_id
			) ;
			// Set table name, if specified
			if( file_spec.size() == 3 ) {
				table_storage->set_table_name( file_spec[2] ) ;
			}
			result = table_storage ;
		} else {
			assert( file_spec[0] == "flat" ) ;
			result = FlatFileOutputter::create(
				file_spec[1],
				analysis_name,
				metadata
			) ;
		}

		return result ;
	}
}

