
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <map>
#include <vector>
#include <string>
#include "genfile/Error.hpp"
#include "statfile/read_values.hpp"

namespace statfile {
	namespace impl {
		std::map< std::size_t, std::size_t > get_indices_of_columns(
 			std::vector< std::string > const& source_column_names,
			std::string const& column_names
		)
			// Get a map from column indices to indices of names in the pipe-delimited string.
			// Names matching the same column can occur more than once.
		{
			std::map< std::size_t, std::size_t > result ;
			std::size_t pos = 0, delim_pos ;
			std::size_t count = 0 ;
			do {
				delim_pos = column_names.find( '|', pos ) ;
				if( delim_pos == std::string::npos ) {
					delim_pos = column_names.size() ;
				}

				std::string column_name = column_names.substr( pos, delim_pos - pos ) ;

				std::vector< std::string >::const_iterator where
					= std::find(
						source_column_names.begin(),
						source_column_names.end(),
						column_name
					) ;

				if( where == source_column_names.end() ) {
					throw genfile::KeyNotFoundError( "source_column_names argument to statfile::impl::get_indices_of_columns()", column_name ) ;
				}

				std::size_t column_index = std::size_t( where - source_column_names.begin() ) ;

				if( result.find( column_index ) != result.end() ) {
					throw genfile::DuplicateKeyError( "source_column_names argument to statfile::impl::get_indices_of_columns()", column_name ) ;
				}
				else {
					result[ column_index ] = count++ ;
				}

				pos = delim_pos + 1 ;
			}
			while( delim_pos != column_names.size() ) ;
			return result ;
		}
	}
}
