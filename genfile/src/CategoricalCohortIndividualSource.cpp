#include <vector>
#include <set>
#include <fstream>
#include <boost/variant/variant.hpp>
#include <boost/bind.hpp>
#include "genfile/CategoricalCohortIndividualSource.hpp"

namespace genfile {
	// Read a sample file in the format described here:
	// http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format_new.html
	CategoricalCohortIndividualSource::CategoricalCohortIndividualSource( std::string const& filename, std::string const& missing_values ):
		FromFileCohortIndividualSource(
			filename,
			string_utils::split_discarding_empty_entries( missing_values, ",", " \t" ),
			&CategoricalCohortIndividualSource::get_entry_from_string
		)
	{
		check_sample_ids() ;
	}

	CategoricalCohortIndividualSource::CategoricalCohortIndividualSource( std::istream& stream, std::string const& missing_values ):
		FromFileCohortIndividualSource(
			stream,
			string_utils::split_discarding_empty_entries( missing_values, ",", " \t" ),
			&CategoricalCohortIndividualSource::get_entry_from_string
		)
	{
		check_sample_ids() ;
	}

	CohortIndividualSource::Entry CategoricalCohortIndividualSource::get_entry_from_string( std::string const& entry_as_string, ColumnType column_type ) {
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
				return Entry( entry_as_string ) ;
				break ;
			}
			case e_BINARY_PHENOTYPE: {
				int value ;
				if( entry_as_string == "control" ) {
					value = 0 ;
				}
				else if( entry_as_string == "case" ) {
					value = 1 ;
				}
				else {
					value = string_utils::to_repr< int >( entry_as_string ) ;
					if( value != 0 && value != 1 ) {
						throw string_utils::StringConversionError() ;
					}
				}
				return Entry( value ) ;
				break ;
			}
			default:
				assert(0) ;
		}
	}
	
	void CategoricalCohortIndividualSource::check_sample_ids() const {
		// combination of id_1 and id_2 fields should be unique.
		{
			std::set< std::string > sample_ids ;
			for( std::size_t i = 0; i < get_number_of_individuals(); ++i ) {
				std::string id = get_entry( i, "id_1" ).as< std::string >() + " " + get_entry( i, "id_2" ).as< std::string >();
				if( sample_ids.find( id ) != sample_ids.end() ) {
					throw DuplicateIndividualError( get_filename(), get_entry( i, "id_1" ).as< std::string >(), get_entry( i, "id_2" ).as< std::string >(), i + 2) ;
				}
				sample_ids.insert( id ) ;
			}
		}
	}
}
