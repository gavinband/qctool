#ifndef RFORMATSTATSOURCE_HPP
#define RFORMATSTATSOURCE_HPP

#include <vector>
#include <string>
#include "statfile/StatSource.hpp"
#include "statfile/IstreamAggregator.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"

namespace statfile {
	// Outputs numerical data in a format suitable for reading with
	// R's read.table().
	class RFormatStatSource: public ColumnNamingStatSource< BuiltInTypeStatSource >, public IstreamAggregator
	{
		typedef ColumnNamingStatSource< BuiltInTypeStatSource > base_t ;
	public:
		RFormatStatSource( std::string const& filename ) ;
		RFormatStatSource( std::auto_ptr< std::istream > stream_ptr ) ;
	public:
		std::size_t number_of_rows() const { return m_number_of_rows ; } ;

		std::string get_descriptive_text() const { return m_descriptive_text ; }

		void reset_to_start() ;

	protected:

		void read_value( int32_t& ) ;
		void read_value( uint32_t& ) ;
		void read_value( std::string& ) ;
		void read_value( double& ) ;
		void ignore_value() ;
		void ignore_all() ;

	private:

		void setup( std::auto_ptr< std::istream > stream_ptr ) ;
		void setup( std::string const& filename ) ;

		std::string read_descriptive_comments() ;
		void read_column_names() ;
		std::size_t count_remaining_lines() ;

	private:
		char m_comment_character ;
		std::size_t m_number_of_comment_lines ;
		std::size_t m_number_of_rows ;
		std::string m_descriptive_text ;
	} ;
}

#endif
