#ifndef GENFILE_FROMFILECOHORTINDIVIDUALSOURCE_HPP
#define GENFILE_FROMFILECOHORTINDIVIDUALSOURCE_HPP

#include <boost/variant.hpp>
#include <boost/function.hpp>
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
		std::size_t get_number_of_covariates() const ;
		std::size_t get_number_of_phenotypes() const ;

	public:
		// Return the column names and types
		std::vector< SingleColumnSpec > get_column_spec() const ;

		// Return the entry in the specified row and column
		Entry get_entry( std::size_t sample_i, std::string const& column_name ) const ;

		// Return the filename from which we read our information, or "(none)" if a stream was supplied. 
		std::string const& get_filename() const ;

	protected:
		std::size_t get_index_of_column( std::string const& column_name ) const ;

	private:
		std::string const m_filename ;
		std::vector< std::string > const m_missing_values ;
		GetEntryFromString m_get_entry_from_string ;
		std::vector< std::string > m_column_names ;
		std::vector< ColumnType > m_column_types ;
		// Entries stored by sample and then by column
		std::vector< std::vector< Entry > > m_entries ;

	protected:
		std::size_t find_column_name( std::string const& column_name ) const ;
	private:
		void setup( std::istream& str ) ;
		void unsafe_setup( std::istream& stream ) ;
		std::vector< std::string > read_column_names( std::istream& stream ) const ;
		std::vector< ColumnType > read_column_types( std::istream& stream, std::vector< std::string > const& column_names ) const ;
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
