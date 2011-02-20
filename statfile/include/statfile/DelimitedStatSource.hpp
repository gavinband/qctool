#ifndef DELIMITEDSTATSOURCE_HPP
#define DELIMITEDSTATSOURCE_HPP

#include <vector>
#include <string>
#include <sstream>
#include "statfile/StatSource.hpp"
#include "statfile/IstreamAggregator.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"

namespace statfile {
	
	class DelimitedStatSource: public ColumnNamingStatSource< BuiltInTypeStatSource >, public IstreamAggregator
		// Class DelimitedStatSource will read files with fields separated by spaces, tabs, commas, or other characters.
		// It allows a number of comment lines at the start of the file, followed by a list of column names (one per column).
		// Although it will strip quotes from around fields, this is currently incorrectly implemented in that
		// the splitting character may not occur within the quoted strings.
	{
		typedef ColumnNamingStatSource< BuiltInTypeStatSource > base_t ;
	public:
		DelimitedStatSource( std::string const& filename, std::string delimiter ) ;
		DelimitedStatSource( std::auto_ptr< std::istream > stream_ptr, std::string delimiter ) ;

		void reset_to_start() ;

	public:
		std::size_t number_of_rows() const { return m_number_of_rows ; } ;

		std::string const& get_descriptive_text() const { return m_descriptive_text ; }

	protected:

		void read_value( int32_t& ) ;
		void read_value( uint32_t& ) ;
		void read_value( std::string& ) ;
		void read_value( double& ) ;
		void ignore_value() ;
		void ignore_all() ;

	private:
		template< typename T >
		void do_read_value( T& value ) {
			if( current_column() == 0 ) {
				read_one_line() ;
			}
			std::istringstream istr( m_current_fields[ current_column() ]) ;
			istr >> value ;
			istr.peek() ;
			if( !istr.eof() ) {
				throw genfile::MalformedInputError( get_source_spec(), number_of_rows_read(), current_column() ) ;
			}
		}

		void read_one_line() ;
		std::vector< std::string > split_line( std::string const& line, std::string const& delimiter, std::string const& strip_chars ) const ;
		static std::string strip( std::string const& string_to_strip, std::string const& strip_chars ) ;
		static std::string get_unquoted_substring( std::string const& big_string, std::size_t pos, std::size_t length, std::string const& quotes ) ;

		void setup( std::auto_ptr< std::istream > stream_ptr ) ;
		void setup( std::string const& filename ) ;

		std::string read_descriptive_comments() ;
		void read_column_names() ;
		std::size_t count_remaining_lines() ;
		std::string get_source_spec() const ;
	private:
		std::string const m_filename ;
		std::vector< std::string > m_current_fields ;
		char const m_comment_character ;
		std::string const m_delimiter ;
		std::string const m_quotes ;
		std::size_t m_number_of_rows ;
		std::string m_descriptive_text ;
		std::size_t m_number_of_comment_lines ;
	} ;
	
	// Specialisation of do_read_value for strings, since we don't need the istringstream here.
	template<>
	void DelimitedStatSource::do_read_value< std::string >( std::string& value ) ;

	// Specialisation of do_read_value for doubles, to deal with infinities
	// represented in the file as "inf".
	template<>
	void DelimitedStatSource::do_read_value< double >( double& value ) ;
}

#endif
