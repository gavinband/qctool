#ifndef STATFILE_UTILS_HPP
#define STATFILE_UTILS_HPP

#include <vector>
#include <map>
#include <string>

namespace statfile {
	
	enum CompressionType { e_NoCompression = 0, e_GzipCompression = 1 } ;
	enum FileFormatType { e_UnknownFormat = 0, e_RFormat = 1 } ;

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

	struct IgnoreOne {
	} ;
	IgnoreOne ignore() ;

	struct IgnoreAll {
	} ;
	IgnoreAll ignore_all() ;

	struct EndRow {	
	} ;
	EndRow end_row() ;

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
	CompressionType get_compression_type_indicated_by_filename( std::string const& filename ) ;

	std::string strip_file_format_extension_if_present( std::string const& filename ) ;

	std::auto_ptr< std::istream > open_text_file_for_input( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::ostream > open_text_file_for_output( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::istream > open_binary_file_for_input( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::ostream > open_binary_file_for_output( std::string filename, CompressionType compression_type ) ;

	struct StatError: public std::exception { char const* what() const throw() { return "StatError" ; } } ;
	struct FileNotOpenedError: public StatError { char const* what() const throw() { return "FileNotOpenedError" ; } } ;
	struct FormatUnsupportedError: public StatError { char const* what() const throw() { return "FormatUnsupportedError" ; } } ;
	struct FileStructureInvalidError: public StatError { char const* what() const throw() { return "FileStructureInvalidError" ; } } ;
	struct FileHasTwoTrailingNewlinesError: public StatError { char const* what() const throw() { return "FileHasTwoTrailingNewlinesError" ; } } ;
}

#endif