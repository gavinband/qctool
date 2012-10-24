
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef STATFILE_DelimitedStatSink_STATSINK_HPP
#define STATFILE_DelimitedStatSink_STATSINK_HPP

#include <vector>
#include <string>
#include <memory>
#include "genfile/Chromosome.hpp"
#include "genfile/GenomePosition.hpp"
#include "genfile/MissingValue.hpp"
#include "statfile/StatSink.hpp"
#include "statfile/OstreamAggregator.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"

namespace statfile {
	// Outputs numerical data in a format suitable for reading with
	// R's read.table().
	class DelimitedStatSink: public ColumnNamingStatSink< BuiltInTypeStatSink >, public OstreamAggregator
	{
	public:
		typedef std::auto_ptr< DelimitedStatSink > UniquePtr ;

		DelimitedStatSink( std::string const& filename, std::string const& delimiter ) ;
		DelimitedStatSink( std::auto_ptr< std::ostream > stream_ptr, std::string const& delimiter ) ;

		operator void*() const { return OstreamAggregator::operator void*() ; }

		void set_descriptive_text( std::string const& ) ;
		void write_metadata( std::string const& metadata ) { set_descriptive_text( metadata ) ; }

	protected:

		void write_value( int32_t const& value ) {
			write_value_impl< int32_t >( value ) ;
		}
		void write_value( uint32_t const& value ) {
			write_value_impl< uint32_t >( value ) ;
		}
		void write_value( int64_t const& value ) {
			write_value_impl< int64_t >( value ) ;
		}
		void write_value( uint64_t const& value ) {
			write_value_impl< uint64_t >( value ) ;
		}
		void write_value( std::string const& value ) {
			write_value_impl< std::string >( value ) ;
		}
		void write_value( genfile::GenomePosition const& value ) {
			write_value_impl< genfile::GenomePosition >( value ) ;
		}
		void write_value( genfile::Chromosome const& value ) {
			write_value_impl< genfile::Chromosome >( value ) ;
		}
		void write_value( genfile::MissingValue const& value ) {
			write_value_impl< genfile::MissingValue >( value ) ;
		}
		void write_value( double const& ) ;

	private:
		
		template< typename T >
		void write_value_impl( T const& value ) {
			write_header_if_necessary() ;
			write_seperator_if_necessary() ;
			stream() << value  ;
		}

	private:
	
		void setup( std::auto_ptr< std::ostream > stream_ptr ) ;
		void setup( std::string const& filename ) ;

		void write_header_if_necessary() {
			if( !(state() & e_HaveWrittenSomeData) ) {
				write_descriptive_text() ;
				write_column_names() ;
				stream() << '\n' ;
			}
		}

		void write_seperator_if_necessary() {
			if( current_column() > 0u ) {
				stream() << m_delimiter ;
			}
		}

		void write_descriptive_text() ;
		void write_column_names() ;

		void move_to_next_row_impl() {
			stream() << "\n" ;
		}

	private:
		char m_comment_character ;
		std::string const m_delimiter ;
		bool have_written_header ;
		int m_precision ;
		std::string m_descriptive_text ;
	} ;
}

#endif
