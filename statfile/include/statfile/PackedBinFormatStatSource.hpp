#ifndef PACKEDBINFORMATSTATSOURCE_HPP
#define PACKEDBINFORMATSTATSOURCE_HPP

#include <vector>
#include <string>
#include "statfile/StatSource.hpp"
#include "statfile/IstreamAggregator.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"

namespace statfile {
	// Outputs numerical data in a format suitable for reading with
	// R's read.table().
	class PackedBinFormatStatSource: public ColumnNamingStatSource< BuiltInTypeStatSource >, public IstreamAggregator
	{
		struct FormatInvalidError: public StatError { char const* what() const throw() { return "PackedBinFormatStatSource::FormatInvalidError" ; } } ;
		
		typedef ColumnNamingStatSource< BuiltInTypeStatSource > base_t ;
	public:
		PackedBinFormatStatSource( std::string const& filename ) ;
		PackedBinFormatStatSource( std::auto_ptr< std::istream > stream_ptr ) ;

	public:
		std::size_t number_of_rows() const ;
		std::string const& get_descriptive_text() const ;

		void reset_to_start() ;

	protected:

		void read_value( int32_t& ) ;
		void read_value( uint32_t& ) ;
		void read_value( std::string& ) ;
		void read_value( double& ) ;
		void ignore_value() ;
		void ignore_all() ;

	private:

		void get_next_entry_location() ;
		void read_string_value( std::string& ) ;
		void ignore_next_value() ;

		void setup( std::auto_ptr< std::istream > stream_ptr ) ;
		void setup( std::string const& filename ) ;

		void read_header_block() ;

		// Read a type specifier, without removing it from the stream.
		char read_type() ;
		// Read a type specifier, throw if it is not the expected one, otherwise remove it from stream.
		void read_type( char expected ) ;
		
	private:
		std::string m_descriptive_text, m_version_text ;
		uint32_t m_number_of_rows ;
		uint32_t m_header_block_length ; // used so we can reset the source to the start.
		uint32_t m_next_entry_rows_read, m_next_entry_column ;
	} ;
}

#endif
