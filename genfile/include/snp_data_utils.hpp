#ifndef SNP_DATA_UTILS_HPP
#define SNP_DATA_UTILS_HPP

#define GENFILE_USE_FAST_PARSE_METHODS 1

#include <iostream>
#include <cassert>
#include <string>
#include <memory>

namespace genfile {
	
	enum ChromosomeEnum {
		Chromosome1 = 1,
		Chromosome2,
		Chromosome3,
		Chromosome4,
		Chromosome5,
		Chromosome6,
		Chromosome7,
		Chromosome8,
		Chromosome9,
		Chromosome10,
		Chromosome11,
		Chromosome12,
		Chromosome13,
		Chromosome14,
		Chromosome15,
		Chromosome16,
		Chromosome17,
		Chromosome18,
		Chromosome19,
		Chromosome20,
		Chromosome21,
		Chromosome22,
		XYPseudoAutosomalDNA = 253,
		MitochondrialDNA = 254,
		UnidentifiedChromosome = 255,
	} ;
	
	struct Chromosome
	{
		Chromosome()
			: m_chromosome_e( UnidentifiedChromosome )
		{}

		Chromosome( ChromosomeEnum chromosome_e )
			: m_chromosome_e( chromosome_e )
		{}

		Chromosome( unsigned char c )
			: m_chromosome_e( ChromosomeEnum( c ) )
		{}

		Chromosome( std::string const& chromosome_str ) ;
		
		Chromosome& operator=( ChromosomeEnum chromosome_e ) {
			m_chromosome_e = chromosome_e ;
			return *this ;
		}
		
		bool operator==( Chromosome other ) {
			return m_chromosome_e == other.m_chromosome_e ;
		}
		
		bool operator<=( Chromosome other ) {
			return m_chromosome_e <= other.m_chromosome_e ;
		}

		bool operator>=( Chromosome other ) {
			return m_chromosome_e >= other.m_chromosome_e ;
		}

		bool operator<( Chromosome other ) {
			return m_chromosome_e < other.m_chromosome_e ;
		}

		bool operator>( Chromosome other ) {
			return m_chromosome_e > other.m_chromosome_e ;
		}
		
		operator ChromosomeEnum() const { return m_chromosome_e ; }

	private:
		ChromosomeEnum m_chromosome_e ;
	} ;
	
	std::ostream& operator<<( std::ostream&, Chromosome const& ) ;
	std::istream& operator>>( std::istream&, Chromosome& ) ;
	
	struct SNPDataBase {} ;
	
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

	struct SNPDataError: public std::exception { char const* what() const throw() { return "SNPDataError" ; } } ;
	struct FileNotOpenedError: public SNPDataError { char const* what() const throw() { return "FileNotOpenedError" ; } } ;
	struct FormatUnsupportedError: public SNPDataError { char const* what() const throw() { return "FormatUnsupportedError" ; } } ;
	struct FileStructureInvalidError: public SNPDataError { char const* what() const throw() { return "FileStructureInvalidError" ; } } ;
	struct FileHasTwoConsecutiveNewlinesError: public SNPDataError { char const* what() const throw() { return "FileHasTwoConsecutiveNewlinesError" ; } } ;
	struct ChromosomeNotRecognisedError: public SNPDataError
	{
		ChromosomeNotRecognisedError( std::string const& input )
			: m_input( input )
		{}
		
		virtual ~ChromosomeNotRecognisedError() throw() {} ;
		
		char const* what() const throw() { return "ChromosomeNotRecognisedError" ; }
		
		std::string const& input() const { return m_input ; }
	private:
		
		std::string m_input ;
	} ;
	
	struct ChromosomeMismatchError: public SNPDataError
	{
		ChromosomeMismatchError( Chromosome expected, Chromosome got )
			: m_expected( expected ),
			m_got( got )
		{}
		
		virtual ~ChromosomeMismatchError() throw() {}
		char const* what() const throw() { return "ChromosomeMismatchError" ; }

		Chromosome const& expected() const { return m_expected ; }
		Chromosome const& got() const { return m_got ; }
	private:
		Chromosome m_expected, m_got ;
	} ;

}

#endif