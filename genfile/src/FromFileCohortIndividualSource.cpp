#include <vector>
#include <set>
#include <fstream>
#include <boost/variant/variant.hpp>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/FromFileCohortIndividualSource.hpp"
#include "string_utils/string_utils.hpp"

namespace genfile {
	// Read a sample file in the format described here:
	// http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format_new.html
	FromFileCohortIndividualSource::FromFileCohortIndividualSource(
		std::string const& filename,
		std::vector< std::string > const& missing_values,
		GetEntryFromString get_entry_from_string,
		GetColumnTypeFromString get_column_type_from_string
	):
		m_filename( filename ),
		m_missing_values( missing_values ),
		m_get_entry_from_string( get_entry_from_string ),
		m_get_column_type_from_string( get_column_type_from_string )
	{
		assert( m_get_entry_from_string ) ;
		std::ifstream stream( m_filename.c_str() ) ;
		// stream.exceptions( std::ios::failbit | std::ios::badbit ) ;
		setup( stream ) ; 
	}

	FromFileCohortIndividualSource::FromFileCohortIndividualSource(
		std::istream& stream,
		std::vector< std::string > const& missing_values,
		GetEntryFromString get_entry_from_string,
		GetColumnTypeFromString get_column_type_from_string
	):
		m_filename( "(none)" ),
		m_missing_values( missing_values ),
		m_get_entry_from_string( get_entry_from_string ),
		m_get_column_type_from_string( get_column_type_from_string )
	{
		assert( m_get_entry_from_string ) ;
		setup( stream ) ; 
	}

	std::size_t FromFileCohortIndividualSource::get_number_of_individuals() const { return m_entries.size() ; }

	std::size_t FromFileCohortIndividualSource::get_number_of_covariates() const {
		return std::count( m_column_types.begin(), m_column_types.end(), e_DISCRETE_COVARIATE )
			+ std::count( m_column_types.begin(), m_column_types.end(), e_CONTINUOUS_COVARIATE ) ;
	}

	std::size_t FromFileCohortIndividualSource::get_number_of_phenotypes() const {
		return std::count( m_column_types.begin(), m_column_types.end(), e_BINARY_PHENOTYPE )
			+ std::count( m_column_types.begin(), m_column_types.end(), e_CONTINUOUS_PHENOTYPE ) ;
	}

	FromFileCohortIndividualSource::Entry FromFileCohortIndividualSource::get_entry( std::size_t sample_i, std::string const& column_name ) const {
		std::size_t const column_i = find_column_name( column_name ) ;
		return m_entries[ sample_i ][ column_i ] ;
	}

	std::string const& FromFileCohortIndividualSource::get_filename() const { return m_filename ; }

	std::vector< CohortIndividualSource::SingleColumnSpec > FromFileCohortIndividualSource::get_column_spec() const {
		assert( m_column_names.size() == m_column_types.size() ) ;
		std::vector< SingleColumnSpec > result( m_column_names.size() ) ;
		for( std::size_t i = 0; i < result.size(); ++i ) {
			result[i] = SingleColumnSpec( m_column_names[i], m_column_types[i] ) ;
		}
		return result ;
	}

	std::size_t FromFileCohortIndividualSource::find_column_name( std::string const& column_name ) const {
		std::vector< std::string >::const_iterator where = std::find( m_column_names.begin(), m_column_names.end(), column_name ) ;
		assert( where != m_column_names.end() ) ;
		return std::size_t( where - m_column_names.begin() ) ;
	}

	void FromFileCohortIndividualSource::setup( std::istream& str ) {
		//try {
			unsafe_setup( str ) ;
		//}
		//catch( std::ios::failure const& e ) {
		//	throw InputError( m_filename ) ;
		//}
	}

	void FromFileCohortIndividualSource::unsafe_setup( std::istream& stream ) {
		m_column_names = read_column_names( stream ) ;
		m_column_types = read_column_types( stream, m_column_names ) ;
		assert( m_column_names.size() == m_column_types.size() ) ;
		m_entries = read_entries( stream, m_column_types) ;
		assert( stream.eof() ) ;
	}
	
	std::vector< std::string > FromFileCohortIndividualSource::read_column_names( std::istream& stream ) const {
		std::vector< std::string > result ;
		std::string line ;
		std::getline( stream, line ) ;
		result = string_utils::split_discarding_empty_entries( string_utils::to_lower( line ), " " ) ;
		if( result.size() < 3 ) {
			throw MalformedInputError( m_filename, 0 ) ;
		}
		if( result[0] != "id_1" ) {
			throw MalformedInputError( m_filename, 0, 0 ) ;
		}
		if( result[1] != "id_2" ) {
			throw MalformedInputError( m_filename, 0, 1 ) ;
		}
		if( result[2] != "missing" ) {
			throw MalformedInputError( m_filename, 0, 2 ) ;
		}
		return result ;
	}

	std::vector< CohortIndividualSource::ColumnType > FromFileCohortIndividualSource::read_column_types( std::istream& stream, std::vector< std::string > const& column_names ) const {
		std::vector< ColumnType > result ;
		std::string line ;
		std::getline( stream, line ) ;
		std::vector< std::string > type_strings = string_utils::split_discarding_empty_entries( line, " " ) ;
		if( type_strings.size() != column_names.size() ) {
			throw MalformedInputError( m_filename, 1 ) ;
		}
		for( std::size_t i = 0; i < type_strings.size(); ++i ) {
			try {
				result.push_back( m_get_column_type_from_string( type_strings[i] )) ;
			}
			catch( MalformedInputError const& e ) {
				throw MalformedInputError( m_filename, 1, i ) ;
			}
		}
		assert( result.size() == column_names.size() ) ;

		// check that only the first three columns have type '0'
		for( std::size_t i = 0; i < result.size(); ++i ) {
			if( i < 3 ) {
				if( result[i] != e_ID_COLUMN ) {
					throw MalformedInputError( m_filename, 1 ) ;
				}
			}
			else if( result[i] == e_ID_COLUMN ) {
				throw MalformedInputError( m_filename, 1 ) ;
			}
		}
		result[2] = e_MISSINGNESS_COLUMN ;
		
		return result ;
	}

	CohortIndividualSource::ColumnType FromFileCohortIndividualSource::get_column_type_from_string_strict( std::string const& a_string ) {
		if( a_string.size() != 1 ) {
			throw MalformedInputError() ;
		}
		switch( a_string[0] ) {
			case '0':
				return e_ID_COLUMN ;
				break ;
			case 'D':
				return e_DISCRETE_COVARIATE ;
				break ;
			case 'C':
				return e_CONTINUOUS_COVARIATE ;
				break ;
			case 'B':
				return e_BINARY_PHENOTYPE ;
				break ;
			case 'P':
				return e_CONTINUOUS_PHENOTYPE ;
				break ;
			default:
				throw MalformedInputError() ;
				break ;
		}
	}

	CohortIndividualSource::ColumnType FromFileCohortIndividualSource::get_column_type_from_string_relaxed( std::string const& a_string ) {
		if( a_string.size() != 1 ) {
			throw MalformedInputError() ;
		}
		switch( a_string[0] ) {
			case '0':
				return e_ID_COLUMN ;
				break ;
			case 'D':
			case '1':
			case '2':
				return e_DISCRETE_COVARIATE ;
				break ;
			case 'C':
			case '3':
				return e_CONTINUOUS_COVARIATE ;
				break ;
			case 'B':
				return e_BINARY_PHENOTYPE ;
				break ;
			case 'P':
				return e_CONTINUOUS_PHENOTYPE ;
				break ;
			default:
				throw MalformedInputError() ;
				break ;
		}
	}
	

	std::vector< std::vector< CohortIndividualSource::Entry > > FromFileCohortIndividualSource::read_entries( std::istream& stream, std::vector< ColumnType > const& column_types ) const {
		std::vector< std::vector< Entry > > result ;
		std::string line ;
		while( std::getline( stream, line ) ) {
			std::vector< std::string > string_entries = string_utils::split_discarding_empty_entries( line, " " ) ;
			if( string_entries.size() != column_types.size() ) {
				throw MalformedInputError( m_filename, 2 + result.size() ) ;
			}
			result.push_back( get_checked_entries( string_entries, column_types, 2 + result.size() )) ;
		}
		return result ;
	}

	std::vector< CohortIndividualSource::Entry > FromFileCohortIndividualSource::get_checked_entries(
		std::vector< std::string > const& string_entries,
		std::vector< ColumnType > const& column_types,
		std::size_t line_number
	) const {
		std::vector< Entry > result ;
		assert( string_entries.size() == column_types.size() ) ;
		for( std::size_t i = 0 ; i < string_entries.size(); ++i ) {
			try {
				result.push_back( get_possibly_missing_entry_from_string( string_entries[i], column_types[i] )) ;
			}
			catch( string_utils::StringConversionError const& e ) {
				throw MalformedInputError( m_filename, line_number, i ) ;
			}
			catch( UnexpectedMissingValueError const& e ) {
				throw UnexpectedMissingValueError( m_filename, line_number, i ) ;
			}
		}
		return result ;
	}

	CohortIndividualSource::Entry FromFileCohortIndividualSource::get_possibly_missing_entry_from_string( std::string const& entry_as_string, ColumnType column_type ) const {
		if( std::find( m_missing_values.begin(), m_missing_values.end(), entry_as_string ) != m_missing_values.end() ) {
			if( column_type == e_ID_COLUMN ) {
				// Missing values are not allowed in ID columns.
				throw UnexpectedMissingValueError( "(none)", 0, 0 ) ;
			}
			else {
				return Entry( MissingValue() ) ;
			}
		}
		else {
			return m_get_entry_from_string( entry_as_string, column_type ) ;
		}
	}
}
