
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef BINFORMATSTATSINK_HPP
#define BINFORMATSTATSINK_HPP

#include <vector>
#include <string>
#include "statfile/StatSink.hpp"
#include "statfile/OstreamAggregator.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"

namespace statfile {
	// Outputs numerical data in a format suitable for reading with
	// R's read.table().
	class BinFormatStatSink: public ColumnNamingStatSink< BuiltInTypeStatSink >, public OstreamAggregator
	{
		typedef ColumnNamingStatSink< BuiltInTypeStatSink > base_t ;
	public:
		// Construct from a filename
		BinFormatStatSink( std::string const& filename, uint32_t number_of_id_columns = 3 ) ;
		// Construct from a stream (must be opened in binary mode).
		BinFormatStatSink( std::auto_ptr< std::ostream > stream_ptr, uint32_t number_of_id_columns = 3 ) ;
		// Finalise written file
		~BinFormatStatSink() ;

		operator bool() const { return OstreamAggregator::operator bool() ; }

		void set_descriptive_text( std::string const& ) ;

	protected:

		void write_value( int32_t const& value ) ;
		void write_value( uint32_t const& value ) ;
		void write_value( std::string const& value ) ;
		void write_value( double const& ) ;

	private:

		void add_column_impl( std::string const& name ) ;

		void setup( std::auto_ptr< std::ostream > stream_ptr ) ;
		void setup( std::string const& filename ) ;

		uint32_t get_header_block_length() const ;
		void write_header_block() ;

	private:
		uint32_t m_number_of_id_columns ;
		std::string m_descriptive_text, m_version_text ;
	} ;
}

#endif
