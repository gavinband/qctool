
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef STATFILE_UTILS_HPP
#define STATFILE_UTILS_HPP

#include <vector>
#include <map>
#include <string>

#include "genfile/FileUtils.hpp"

namespace statfile {
	
	enum FileFormatType {
		e_UnknownFormat = 0,
		e_SpaceDelimited = 1,					// Readable by R's read.table( <filename>, header = T )
		e_CommaDelimitedFormat = 2,		// CSV format, same number of column headers as columns.
		e_TabDelimitedFormat = 3		// Tab-delimited format, same number of column headers as columns.
	} ;

	struct MapValueSetter
	{
		MapValueSetter( std::map< std::string, std::string >& map ) ;
		void operator()( std::string const& key, std::string const& value ) ;
	private:
		std::map< std::string, std::string > m_map ;
	} ;

	template< typename T >
	struct ValueSetter
	{
		ValueSetter( T& t ): m_t(t) {}
		void operator()( T const& t ) { m_t = t ; }
	private:
		T& m_t ;
	} ;
	
	template< typename T > ValueSetter< T > set( T& t ) {
		return ValueSetter< T >( t ) ;
	}

	struct IgnoreSome {
		IgnoreSome( std::size_t n = 1 ): m_number_to_ignore( n ) {}
		std::size_t number_to_ignore() const { return m_number_to_ignore ; }
	private:
		std::size_t m_number_to_ignore ;
	} ;

	IgnoreSome ignore( std::size_t n = 1 ) ;

	struct IgnoreAll {
	} ;
	IgnoreAll ignore_all() ;

	struct EndRow {	
	} ;
	EndRow end_row() ;

	struct BeginData {} ;
	BeginData begin_data() ;
	
	struct ColumnSpec
	{
		ColumnSpec( std::string name, std::size_t repeats ) ;
 
		std::string const& name() const { return m_name ; }
		std::size_t number_of_repeats() const { return m_number_of_repeats ; }

	private:
		std::string m_name ;
		std::size_t m_number_of_repeats ;
	} ;

	typedef std::vector< ColumnSpec > ColumnSpecList ;

	ColumnSpecList operator|( ColumnSpec const& left, ColumnSpec const& right ) ;
	ColumnSpecList operator|( ColumnSpecList const& left, ColumnSpec const& right ) ;

	FileFormatType get_file_format_type_indicated_by_filename( std::string const& filename ) ;

	using genfile::get_compression_type_indicated_by_filename ;
	using genfile::open_text_file_for_input ;
	using genfile::open_text_file_for_input ;
	using genfile::open_text_file_for_output  ;
	using genfile::open_binary_file_for_input ;
	using genfile::open_binary_file_for_output ;

	struct StatError: public std::exception { char const* what() const throw() { return "StatError" ; } } ;
	struct FileNotOpenedError: public StatError {
		FileNotOpenedError( std::string const& filename = "(unknown)" ): m_filename( filename ) {}
		~FileNotOpenedError() throw() {}
		char const* what() const throw() { return "FileNotOpenedError" ; }
		std::string const& filename() const { return m_filename ; }
	private:
		std::string const m_filename ;
	} ;
	struct FormatUnsupportedError: public StatError { char const* what() const throw() { return "FormatUnsupportedError" ; } } ;
	struct FileStructureInvalidError: public StatError { char const* what() const throw() { return "FileStructureInvalidError" ; } } ;
	struct FileHasTwoConsecutiveNewlinesError: public StatError { char const* what() const throw() { return "FileHasTwoConsecutiveNewlinesError" ; } } ;
	struct InvalidColumnNameError: public StatError
	{
		InvalidColumnNameError( std::string const& name ): m_name( name ) {}
		~InvalidColumnNameError() throw() {}
		char const* what() const throw() { return "InvalidColumnNameError" ; }
		std::string const& name() const { return m_name ; }
	private:
		std::string const m_name ;
	} ;
}

#endif
