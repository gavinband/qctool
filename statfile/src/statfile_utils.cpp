
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <map>
#include <string>
#include <cstdio>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include "genfile/FileUtils.hpp"
#include "genfile/Error.hpp"
#include "statfile/statfile_utils.hpp"

namespace statfile {
	IgnoreSome ignore( std::size_t n ) { return IgnoreSome( n ) ; }
	IgnoreAll ignore_all() { return IgnoreAll() ; }
	EndRow end_row() { return EndRow() ; }
	RestartRow restart_row() { return RestartRow() ; }
	BeginData begin_data() { return BeginData() ; }

	ColumnSpec::ColumnSpec( std::string name, std::size_t repeats ) 
		: m_name( name ), m_number_of_repeats( repeats )
	{
		assert( m_number_of_repeats > 0 ) ;
	}

	std::vector< ColumnSpec > operator|( ColumnSpec const& left, ColumnSpec const& right ) {
		std::vector< ColumnSpec > result ;
		result.reserve(2) ;
		result.push_back( left ) ;
		result.push_back( right ) ;
		return result ;
	}

	std::vector< ColumnSpec > operator|( std::vector< ColumnSpec > const& left, ColumnSpec const& right ) {
		std::vector< ColumnSpec > result = left ;
		result.push_back( right ) ;
		return result ;
	}
	
	FileFormatType get_file_format_type_indicated_by_filename( std::string const& filename ) {
		std::map< std::string, FileFormatType > types  ;
		types[ ".txt" ]     = types[ ".txt.gz" ]        = e_SpaceDelimited ;
		types[ ".csv" ]     = types[ ".csv.gz" ] 		= e_CommaDelimitedFormat ;
		types[ ".tsv" ]     = types[ ".tsv.gz" ] 		= e_TabDelimitedFormat ;

		for(
			std::map< std::string, FileFormatType >::const_iterator i = types.begin();
			i != types.end();
			++i
		) {
			if(
				filename.size() >= i->first.size() &&
				filename.substr( filename.size()- i->first.size(), i->first.size() ) == i->first
			) {
				return i->second ;
			}
		}
		return e_UnknownFormat ;
	}

	std::string strip_file_format_extension_if_present( std::string const& filename ) {
		std::string recognised_extensions[1] = {
			".rdata"
		} ;

		for( std::size_t i = 0; i < 1u; ++i ) {
			if( filename.substr( filename.size() - recognised_extensions[i].size(), recognised_extensions[i].size() ) == recognised_extensions[i] ) {
				return filename.substr( 0, filename.size() - recognised_extensions[i].size() ) ;
			}
		}

		return filename ;
	}
}
