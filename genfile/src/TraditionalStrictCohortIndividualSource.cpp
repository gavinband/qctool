#include <vector>
#include <set>
#include <fstream>
#include <cassert>
#include <boost/variant/variant.hpp>
#include <boost/bind.hpp>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/FromFileCohortIndividualSource.hpp"
#include "genfile/TraditionalStrictCohortIndividualSource.hpp"

namespace genfile {
	// Read a sample file in the format described here:
	// http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format_new.html
	TraditionalStrictCohortIndividualSource::TraditionalStrictCohortIndividualSource( std::string const& filename, std::string const& missing_values ):
		FromFileCohortIndividualSource(
			filename,
			string_utils::split_discarding_empty_entries( missing_values, ",", " \t" ),
			&TraditionalStrictCohortIndividualSource::get_entry_from_string,
			&FromFileCohortIndividualSource::get_column_type_from_string_strict
		)
	{
		check_sample_ids() ;
	}

	TraditionalStrictCohortIndividualSource::TraditionalStrictCohortIndividualSource( std::istream& stream, std::string const& missing_values ):
		FromFileCohortIndividualSource(
			stream,
			string_utils::split_discarding_empty_entries( missing_values, ",", " \t" ),
			&TraditionalStrictCohortIndividualSource::get_entry_from_string,
			&FromFileCohortIndividualSource::get_column_type_from_string_strict
		)
	{
		check_sample_ids() ;
	}

	CohortIndividualSource::Entry TraditionalStrictCohortIndividualSource::get_entry_from_string( std::string const& entry_as_string, ColumnType column_type ) {
		switch( column_type ) {
			case e_ID_COLUMN:
				return Entry( entry_as_string ) ;
				break ;
			case e_MISSINGNESS_COLUMN:
			case e_CONTINUOUS_COVARIATE:
			case e_CONTINUOUS_PHENOTYPE:
				return Entry( string_utils::to_repr< double >( entry_as_string )) ;
				break ;
			case e_DISCRETE_COVARIATE: {
				int value = string_utils::to_repr< int >( entry_as_string ) ;
				if( value <= 0 ) {
					throw string_utils::StringConversionError() ;
				}
				return Entry( value ) ;
				break ;
			}
			case e_BINARY_PHENOTYPE: {
				int value = string_utils::to_repr< int >( entry_as_string ) ;
				if( value != 0 && value != 1 ) {
					throw string_utils::StringConversionError() ;
				}
				return Entry( value ) ;
				break ;
			}
			default:
				assert(0) ;
		}
	}
	
	void TraditionalStrictCohortIndividualSource::check_sample_ids() const {
		// each of id_1 and id_2 fields should be unique.
		{
			std::set< std::string > sample_ids ;
			for( std::size_t i = 0; i < get_number_of_individuals(); ++i ) {
				std::string id = get_entry( i, "id_1" ).as< std::string >() ;
				if( sample_ids.find( id ) != sample_ids.end() ) {
					throw DuplicateKeyError( get_filename(), i + 2, find_column_name( "id_1" )) ;
				}
				sample_ids.insert( id ) ;
			}
		}

		{
			std::set< std::string > sample_ids ;
			for( std::size_t i = 0; i < get_number_of_individuals(); ++i ) {
				std::string id = get_entry( i, "id_2" ).as< std::string >() ;
				if( sample_ids.find( id ) != sample_ids.end() ) {
					throw DuplicateKeyError( get_filename(), i + 2, find_column_name( "id_2" ) ) ;
				}
				sample_ids.insert( id ) ;
			}
		}
	}
	
}
