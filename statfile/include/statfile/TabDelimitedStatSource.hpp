#ifndef TABDELIMITEDSTATSOURCE_HPP
#define TABDELIMITEDSTATSOURCE_HPP

#include <vector>
#include <string>
#include <sstream>
#include "statfile/StatSource.hpp"
#include "statfile/IstreamAggregator.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"

namespace statfile {
	// Outputs numerical data in a format suitable for reading with
	// R's read.table().
	class TabDelimitedStatSource: public ColumnNamingStatSource< BuiltInTypeStatSource >, public IstreamAggregator
	{
		typedef ColumnNamingStatSource< BuiltInTypeStatSource > base_t ;
	public:
		TabDelimitedStatSource( std::string const& filename ) ;
		TabDelimitedStatSource( std::auto_ptr< std::istream > stream_ptr ) ;

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
				throw FileIsInvalidError() ;
			}
		}

		void read_one_line() ;
		std::vector< std::string > split_line( std::string const& line, char delimiter ) const ;

		void setup( std::auto_ptr< std::istream > stream_ptr ) ;
		void setup( std::string const& filename ) ;

		std::string read_descriptive_comments() ;
		void read_column_names() ;
		std::size_t count_remaining_lines() ;

	private:
		std::vector< std::string > m_current_fields ;
		char m_comment_character ;
		char m_delimiter ;
		std::size_t m_number_of_rows ;
		std::string m_descriptive_text ;
		std::size_t m_number_of_comment_lines ;
	} ;
	
	// Specialisation of do_read_value for strings, since we don't need the istringstream here.
	template<>
	void TabDelimitedStatSource::do_read_value< std::string >( std::string& value ) ;

	// Specialisation of do_read_value for doubles, to deal with infinities
	// represented in the file as "inf".
	template<>
	void TabDelimitedStatSource::do_read_value< double >( double& value ) ;
}

#endif
