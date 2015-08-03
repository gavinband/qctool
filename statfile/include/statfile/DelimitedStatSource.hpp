
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef DELIMITEDSTATSOURCE_HPP
#define DELIMITEDSTATSOURCE_HPP

#include <vector>
#include <string>
#include <sstream>
#include <boost/optional.hpp>
#include "genfile/string_utils/slice.hpp"
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
		typedef ColumnNamingStatSource< BuiltInTypeStatSource > Base ;
		typedef boost::optional< std::string > OptionalString ;
	public:
		DelimitedStatSource(
			std::string const& filename,
			std::string delimiter,
			OptionalString const& ignore_until = OptionalString(),
			OptionalString const& ignore_from = OptionalString()
		) ;
		DelimitedStatSource(
			std::auto_ptr< std::istream > stream_ptr,
			std::string delimiter,
			OptionalString const& ignore_until = OptionalString(),
			OptionalString const& ignore_from = OptionalString()
		) ;

		void reset_to_start() ;
		operator bool() const ;
		
	public:
		OptionalCount number_of_rows() const { return m_number_of_rows ; } ;

		std::string get_descriptive_text() const { return m_descriptive_text ; }

	protected:

		using Base::read_value ;
		void read_value( int32_t& ) ;
		void read_value( uint32_t& ) ;
		void read_value( std::string& ) ;
		void read_value( double& ) ;
		void ignore_value() ;
		void ignore_all() ;
		void move_to_next_row_impl() ;
		void restart_row() ;

	private:
		template< typename T >
		void do_read_value( T& value ) {
			assert( m_have_more_data ) ;
			std::istringstream istr( m_current_fields[ current_column() ]) ;
			istr >> value ;
			istr.peek() ;
			if( !istr.eof() ) {
				throw genfile::MalformedInputError( get_source_spec(), number_of_rows_read(), current_column() ) ;
			}
		}

		void read_one_line() ;
		std::vector< genfile::string_utils::slice > split_line( std::string const& line, std::string const& delimiter, std::string const& strip_chars ) const ;
		static std::string strip( std::string const& string_to_strip, std::string const& strip_chars ) ;
		static genfile::string_utils::slice get_unquoted_substring( std::string const& big_string, std::size_t pos, std::size_t length, std::string const& quotes ) ;

		void setup( std::auto_ptr< std::istream > stream_ptr ) ;
		void setup( std::string const& filename ) ;

		std::string read_comments() ;
		void read_column_names() ;
		std::size_t count_remaining_lines() ;
		std::string get_source_spec() const ;
	private:
		std::string const m_filename ;
		std::string m_current_line ;
		bool m_have_more_data ;
		std::vector< genfile::string_utils::slice > m_current_fields ;
		char const m_comment_character ;
		std::string const m_delimiter ;
		std::string const m_quotes ;
		OptionalCount m_number_of_rows ;
		std::string m_descriptive_text ;
		std::size_t m_number_of_comment_lines ;
		std::size_t m_number_of_ignored_lines ;
		
		// m_ignore_until and m_ignore_from represent a way to make the source
		// ignore a header section and a trailing section of the file.
		// This is useful for the type of csv that Excel sometimes writes.
		boost::optional< std::string > m_ignore_until ;
		boost::optional< std::string > m_ignore_from ;
		
		void reset_stream_to_start() ;
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
