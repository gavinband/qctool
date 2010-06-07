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
	// To allow flexibility, this class takes as a parameters two functions.  The
	// first, get_entry_from_string, should return an Entry object given its string
	// representation.  
	//
	// The second, get_column_type_from_string, returns a ColumnType given a string.  Two possible
	// implementations are provided as static functions; the first implements the file format above
	// while the second also accepts files formatted in "v1" format.
	// It does this by mapping the column types as follows:
	// 1,2 -> discrete covariate, 3 -> continuous covariate, P -> phenotype
	//
	class FromFileCohortIndividualSource: public CohortIndividualSource
	{
	public:
		typedef boost::function< Entry ( std::string const&, ColumnType ) > GetEntryFromString ;
		typedef boost::function< ColumnType ( std::string const& ) > GetColumnTypeFromString ;
		typedef boost::function< void ( CohortIndividualSource const& ) > Check ;
	public:
		FromFileCohortIndividualSource(
			std::string const& filename,
			std::vector< std::string > const& missing_values,
			GetEntryFromString get_entry_from_string,
			GetColumnTypeFromString get_column_type_from_string
		) ;
		FromFileCohortIndividualSource(
			std::istream& stream,
			std::vector< std::string > const& missing_values,
			GetEntryFromString get_entry_from_string,
			GetColumnTypeFromString get_column_type_from_string
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
		GetColumnTypeFromString m_get_column_type_from_string ;
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
	protected:
		// For convenience, here are two possible implementations for the get_column_type_from_string constructor argument
		static ColumnType get_column_type_from_string_strict( std::string const& string ) ;
		static ColumnType get_column_type_from_string_relaxed( std::string const& string ) ;
	} ;
}

#endif
