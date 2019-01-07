
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef BGEN_INDEX_QUERY_HPP
#define BGEN_INDEX_QUERY_HPP

#include <boost/function.hpp>
#include <boost/optional.hpp>
#include <stdint.h>
#include <vector>
#include <string>
#include <ctime>
#include "genfile/db/sqlite3.hpp"
#include "genfile/bgen/Query.hpp"

namespace genfile {
	namespace bgen {
		// Base class representing a query against a BGEN file index
		struct IndexQuery {
		public:
			// We use std::auto_ptr to avoid using C++11 features here.
			typedef std::auto_ptr< IndexQuery > UniquePtr ;
			typedef uint8_t byte_t ;
			typedef boost::function< void ( std::size_t n, boost::optional< std::size_t > total ) > ProgressCallback ;
			struct FileMetadata ;
			typedef boost::optional< FileMetadata > OptionalFileMetadata ;
			typedef Query::GenomicRange GenomicRange ;
			typedef std::pair< int64_t, int64_t> FileRange ;

		public:
			static UniquePtr create(
				std::string const& filename,
				Query const& query,
				ProgressCallback callback = ProgressCallback(),
				std::string const& table_name = "Variant"
			) ;

		public:
			virtual ~IndexQuery() {} ;
			virtual OptionalFileMetadata const& file_metadata() const = 0 ;

		public:
			// Report the number of variants in this query.
			virtual std::size_t number_of_variants() const = 0 ;
			// Report the number of variants in this query.
			virtual FileRange locate_variant( std::size_t index ) const = 0 ;

			struct FileMetadata {
				FileMetadata():
					size(-1)
				{}

				FileMetadata( FileMetadata const& other ):
					filename( other.filename ),
					size( other.size ),
					last_write_time( other.last_write_time ),
					first_bytes( other.first_bytes )
				{}

				FileMetadata& operator=( FileMetadata const& other ) {
					filename = other.filename ;
					size = other.size ;
					last_write_time = other.last_write_time ;
					first_bytes = other.first_bytes ;
					return *this ;
				}

				std::string filename ;
				int64_t size ;
				std::time_t last_write_time ;
				std::vector< byte_t > first_bytes ;
			} ;
		} ;
		
		// Class for index queries implemented using a sqlite file, a la bgenix.
		struct SqliteIndexQuery: public IndexQuery {
		public:
			// We use auto_ptr to avoid using C++11 features here.
			typedef std::auto_ptr< SqliteIndexQuery > UniquePtr ;

		public:
			SqliteIndexQuery( std::string const& filename, Query const& query, ProgressCallback callback = ProgressCallback(), std::string const& table_name = "Variant" ) ;
			OptionalFileMetadata const& file_metadata() const ;
			std::size_t number_of_variants() const ;
			FileRange locate_variant( std::size_t index ) const ;

		private:
			db::Connection::UniquePtr open_connection( std::string const& filename ) const ;
			OptionalFileMetadata load_metadata( db::Connection& connection ) const ;
			db::Connection::StatementPtr build_query() const ;

			// Implement the specified query.
			void implement_query( Query const& query, ProgressCallback callback = ProgressCallback() ) ;
			// Internal query implementation methods
			// Include variants in a range
			SqliteIndexQuery& include_range( GenomicRange const& range ) ;
			// Exclude variants in a range
			SqliteIndexQuery& exclude_range( GenomicRange const& range ) ;
			// Include variants in a range
			SqliteIndexQuery& include_ranges( std::vector< GenomicRange > const& ranges ) ;
			// Exclude variants in a range
			SqliteIndexQuery& exclude_ranges( std::vector< GenomicRange > const& ranges ) ;
			// Include variants with one of the given rsids.  The list provided must be unique.
			SqliteIndexQuery& include_rsids( std::vector< std::string > const& ids ) ;
			// Exclude variants with one of the given rsids.  The list provided must be unique.
			SqliteIndexQuery& exclude_rsids( std::vector< std::string > const& ids ) ;
			// Initialise query with the above
			void initialise( ProgressCallback callback = ProgressCallback() ) ;

		private:
			db::Connection::UniquePtr m_connection ;
			OptionalFileMetadata const m_metadata ;
			std::string const m_index_table_name ;
			struct QueryParts {
				std::string join ;
				std::string inclusion ;
				std::string exclusion ;
			} ;
			QueryParts m_query_parts ;
			bool m_initialised ;
			std::vector< std::pair< int64_t, int64_t> > m_positions ;
		} ;
		
	}
}

#endif
