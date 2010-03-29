#ifndef SNP_DATA_UTILS_HPP
#define SNP_DATA_UTILS_HPP

#define GENFILE_USE_FAST_PARSE_METHODS 1

#include <iostream>
#include <cassert>
#include <string>
#include <memory>

#include "genfile/Chromosome.hpp"

namespace genfile {
	
	struct SNPDataBase {} ;
	
	struct Ignorer
	{
		template< typename T > void operator()( T const& ) {}
		template< typename T > void operator()( T const&, T const& ) {}
		template< typename T > void operator()( T const&, T const&, T const& ) {}
		template< typename T1, typename T2 > void operator()( T1 const&, T2 const&, T2 const&, T2 const& ) {}
	}  ;

	Ignorer ignore() ;

	template< typename T >
	struct ValueSetter
	{
		ValueSetter( T& t ): m_t(t) {}
		template< typename T2 >
		void operator()( T2 const& t ) { m_t = T(t) ; }
	private:
		T& m_t ;
	} ;
	
	template<>
	struct ValueSetter< std::string >
	{
		ValueSetter< std::string >( std::string& t ): m_t(t) {}
		void operator()( char const& c ) { m_t.assign( std::size_t(1), c ) ; }
		template< typename T2 > void operator()( T2 const& t ) { m_t = t ; }
	private:
		std::string& m_t ;
	} ;
/*
	template<>
	struct ValueSetter< Chromosome >
	{
		ValueSetter< Chromosome >( Chromosome& t ): m_t(t) {}
		void operator()( unsigned char c ) { m_t = Chromosome( c ) ; }
		template< typename T2 > void operator()( T2 const& t ) { m_t = t ; }
	private:
		Chromosome& m_t ;
	} ;
*/
	template< typename T > ValueSetter< T > set_value( T& t ) {
		return ValueSetter< T >( t ) ;
	}

	enum CompressionType { e_NoCompression = 0, e_GzipCompression = 1 } ;

	bool filename_indicates_gen_format( std::string const& filename ) ;
	bool filename_indicates_bgen_format( std::string const& filename ) ;
	bool filename_indicates_gen_or_bgen_format( std::string const& filename ) ;
	CompressionType get_compression_type_indicated_by_filename( std::string const& filename ) ;
	Chromosome get_chromosome_indicated_by_filename( std::string const& filename ) ;

	std::string get_gen_file_extension_if_present( std::string const& filename ) ;
	std::string strip_gen_file_extension_if_present( std::string const& filename ) ;

	std::auto_ptr< std::istream > open_text_file_for_input( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::ostream > open_text_file_for_output( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::istream > open_binary_file_for_input( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::ostream > open_binary_file_for_output( std::string filename, CompressionType compression_type ) ;

	std::string create_temporary_filename() ;

	struct SNPDataError: public std::exception { char const* what() const throw() { return "genfile::SNPDataError" ; } } ;
	struct OperationFailedError: public std::exception { char const* what() const throw() { return "genfile::OperationFailedError" ; } } ;
	struct FileNotOpenedError: public SNPDataError { char const* what() const throw() { return "genfile::FileNotOpenedError" ; } } ;
	struct FormatUnsupportedError: public SNPDataError { char const* what() const throw() { return "genfile::FormatUnsupportedError" ; } } ;
	struct FileStructureInvalidError: public SNPDataError { char const* what() const throw() { return "genfile::FileStructureInvalidError" ; } } ;
	struct FileHasTwoConsecutiveNewlinesError: public SNPDataError { char const* what() const throw() { return "genfile::FileHasTwoConsecutiveNewlinesError" ; } } ;
	struct FileNotFoundError: public SNPDataError	
	{
		FileNotFoundError( std::string const& filespec ): m_filespec( filespec ) {}
		~FileNotFoundError() throw() {}
		char const* what() const throw() { return "genfile::FileNotFoundError" ; }
		std::string const& filespec() const { return m_filespec ; }
		private:
			std::string const m_filespec ;
	} ;
}

#endif
