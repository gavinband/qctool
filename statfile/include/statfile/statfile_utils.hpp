#ifndef STATFILE_UTILS_HPP
#define STATFILE_UTILS_HPP

#include <vector>
#include <map>
#include <string>

namespace statfile {
	
	enum CompressionType { e_NoCompression = 0, e_GzipCompression = 1 } ;
	enum FileFormatType {
		e_UnknownFormat = 0,
		e_RFormat = 1,
		e_TabDelimitedFormat = 2,
		e_BinFormat = 3,
		e_PackedBinFormat = 4
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
	std::auto_ptr< std::istream > open_binary_file_for_input( std::string filename, CompressionType compression_type = e_NoCompression ) ;
	std::auto_ptr< std::ostream > open_binary_file_for_output( std::string filename, CompressionType compression_type = e_NoCompression ) ;

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