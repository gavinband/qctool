
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_FROMFILECOHORTINDIVIDUALSOURCE_HPP
#define GENFILE_FROMFILECOHORTINDIVIDUALSOURCE_HPP

#include <map>
#include <string>
#include <boost/variant.hpp>
#include <boost/function.hpp>
#include <boost/optional.hpp>
#include "genfile/CohortIndividualSource.hpp"

namespace genfile {
	// class FromFileCohortIndividualSource
	// This class forms a base class for classes which read sample information from flat files
	// in a format like the one described here:
	// http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format_new.html
	// To allow flexibility, this class takes as a parameters a function, get_entry_from_string,
	// which should return an Entry object given its string representation.
	// This can be used for varying the allowed formats.
	class FromFileCohortIndividualSource: public CohortIndividualSource
	{
	public:
		typedef boost::function< Entry ( std::string const&, ColumnType ) > GetEntryFromString ;
		typedef boost::function< void ( CohortIndividualSource const& ) > Check ;
	public:
		FromFileCohortIndividualSource(
			std::string const& filename,
			std::vector< std::string > const& missing_values,
			GetEntryFromString get_entry_from_string
		) ;
		FromFileCohortIndividualSource(
			std::istream& stream,
			std::vector< std::string > const& missing_values,
			GetEntryFromString get_entry_from_string
		) ;

		virtual ~FromFileCohortIndividualSource() {}

	public:
		std::size_t get_number_of_individuals() const ;

	public:
		// Return the column names and types
		ColumnSpec get_column_spec() const ;

		// Return the entry in the specified row and column
		Entry get_entry( std::size_t sample_i, std::string const& column_name ) const ;

		bool check_for_column( std::string const& column_name ) const ;

		// Return the filename from which we read our information, or "(none)" if a stream was supplied. 
		std::string const& get_filename() const ;
		// Return the filename from which we read our information, or "(unknown)" if a stream was supplied. 
		std::string get_source_spec() const ;

	private:
		std::string const m_filename ;
		std::vector< std::string > const m_missing_values ;
		GetEntryFromString m_get_entry_from_string ;
		std::string m_comments ;
		std::size_t m_number_of_metadata_lines ;
		std::vector< std::string > m_column_names ;
		std::vector< ColumnType > m_column_types ;
		// Entries stored by sample and then by column
		std::vector< std::vector< Entry > > m_entries ;

	protected:
		// Case-insensitive search for a column in the sample file.
		std::size_t find_column_name( std::string const& column_name ) const ;
	private:
		std::vector< std::string >::const_iterator find_column_name_impl( std::string const& column_name ) const ;
		void setup( std::istream& str ) ;
		void unsafe_setup( std::istream& stream ) ;
		std::string read_comments( std::istream& stream ) const ;
		std::vector< std::string > read_column_names( std::istream& stream ) const ;
		boost::optional< std::map< std::string, CohortIndividualSource::ColumnType > > read_column_types_from_comments( std::string const& comments ) const ;
		std::vector< ColumnType > read_column_type_line( std::istream& stream, std::vector< std::string > const& column_names ) const ;
		std::vector< std::vector< Entry > > read_entries( std::istream& stream, std::vector< ColumnType > const& column_types ) const ;
		std::vector< Entry > get_checked_entries(
			std::vector< std::string > const& string_entries,
			std::vector< ColumnType > const& column_types,
			std::size_t line_number
		) const ;
		
		Entry get_possibly_missing_entry_from_string( std::string const& entry_as_string, ColumnType column_type ) const ;
	} ;
}

#endif
