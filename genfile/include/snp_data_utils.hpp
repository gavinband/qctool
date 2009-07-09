#ifndef SNP_DATA_UTILS_HPP
#define SNP_DATA_UTILS_HPP

#define GENFILE_USE_FAST_PARSE_METHODS 1

namespace genfile {
	struct Ignorer
	{
		template< typename T > void operator()( T const& ) {} ;
		template< typename T > void operator()( T const&, T const& ) {} ;
		template< typename T > void operator()( T const&, T const&, T const& ) {} ;
		template< typename T1, typename T2 > void operator()( T1 const&, T2 const&, T2 const&, T2 const& ) {} ;
	} ;

	Ignorer ignore() ;

	template< typename T >
	struct ValueSetter
	{
		ValueSetter( T& t ): m_t(t) {}
		void operator()( T const& t ) { m_t = t ; }
	private:
		T& m_t ;
	} ;

	template< typename T > ValueSetter< T > set_value( T& t ) {
		return ValueSetter< T >( t ) ;
	}
	

	enum CompressionType { e_NoCompression = 0, e_GzipCompression = 1 } ;

	bool filename_indicates_bgen_format( std::string const& filename ) ;
	CompressionType get_compression_type_indicated_by_filename( std::string const& filename ) ;


	std::auto_ptr< std::istream > open_text_file_for_input( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::ostream > open_text_file_for_output( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::istream > open_binary_file_for_input( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::ostream > open_binary_file_for_output( std::string filename, CompressionType compression_type ) ;

	std::string create_temporary_filename() ;

	struct SNPDataError: public std::exception { char const* what() const throw() { return "SNPDataError" ; } } ;
	struct FileNotOpenedError: public SNPDataError { char const* what() const throw() { return "FileNotOpenedError" ; } } ;
	struct FormatUnsupportedError: public SNPDataError { char const* what() const throw() { return "FormatUnsupportedError" ; } } ;
}

#endif