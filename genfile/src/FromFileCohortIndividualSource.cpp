
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <set>
#include <fstream>
#include <boost/variant/variant.hpp>
#include <boost/bind.hpp>
#include "genfile/FileUtils.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/FromFileCohortIndividualSource.hpp"

namespace genfile {
	// Read a sample file in the format described here:
	// http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format_new.html
	FromFileCohortIndividualSource::FromFileCohortIndividualSource(
		std::string const& filename,
		std::vector< std::string > const& missing_values,
		GetEntryFromString get_entry_from_string
	):
		m_filename( filename ),
		m_missing_values( missing_values ),
		m_get_entry_from_string( get_entry_from_string )
	{
		assert( m_get_entry_from_string ) ;
		std::auto_ptr< std::istream > stream = open_text_file_for_input( m_filename ) ;
		// stream.exceptions( std::ios::failbit | std::ios::badbit ) ;
		setup( *stream ) ; 
	}

	FromFileCohortIndividualSource::FromFileCohortIndividualSource(
		std::istream& stream,
		std::vector< std::string > const& missing_values,
		GetEntryFromString get_entry_from_string
	):
		m_filename( "(none)" ),
		m_missing_values( missing_values ),
		m_get_entry_from_string( get_entry_from_string )
	{
		assert( m_get_entry_from_string ) ;
		setup( stream ) ; 
	}

	std::size_t FromFileCohortIndividualSource::get_number_of_individuals() const { return m_entries.size() ; }

	FromFileCohortIndividualSource::Entry FromFileCohortIndividualSource::get_entry( std::size_t sample_i, std::string const& column_name ) const {
		std::size_t const column_i = find_column_name( column_name ) ;
		return m_entries[ sample_i ][ column_i ] ;
	}

	std::string const& FromFileCohortIndividualSource::get_filename() const { return m_filename ; }

	std::string FromFileCohortIndividualSource::get_source_spec() const {
		if( m_filename == "(none)" ) {
			return "(unknown)" ;
		}
		else {
			return "file://" + m_filename ;
		}
	}

	CohortIndividualSource::ColumnSpec FromFileCohortIndividualSource::get_column_spec() const {
		return ColumnSpec( m_column_names, m_column_types ) ;
	}

	bool FromFileCohortIndividualSource::check_for_column( std::string const& column_name ) const {
		std::vector< std::string >::const_iterator where = find_column_name_impl( column_name ) ;
		return ( where != m_column_names.end() ) ;
	}

	std::vector< std::string >::const_iterator FromFileCohortIndividualSource::find_column_name_impl( std::string const& column_name ) const {
		return std::find_if(
			m_column_names.begin(),
			m_column_names.end(),
			boost::bind(
				string_utils::case_insensitive_equality,
				column_name,
				_1
			)
		) ;
	}
	
	std::size_t FromFileCohortIndividualSource::find_column_name( std::string const& column_name ) const {
		std::vector< std::string >::const_iterator where = find_column_name_impl( column_name ) ;
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
		result = string_utils::split_and_strip_discarding_empty_entries( line ) ;
		if( result.size() < 3 ) {
			throw MalformedInputError( m_filename, 0 ) ;
		}
		if( string_utils::to_lower( result[0] ) != "id_1" ) {
			throw MalformedInputError( m_filename, 0, 0 ) ;
		}
		if( string_utils::to_lower( result[1] ) != "id_2" ) {
			throw MalformedInputError( m_filename, 0, 1 ) ;
		}
		if( string_utils::to_lower( result[2] ) != "missing" ) {
			throw MalformedInputError( m_filename, 0, 2 ) ;
		}
		// check for uniqueness
		for( std::size_t i = 0; i < result.size(); ++i ) {
			if( ( std::find( result.begin(), result.end(), result[i] ) - result.begin() ) < i ) {
				throw DuplicateKeyError( m_filename, result[i] ) ;
			}
		}
		return result ;
	}

	std::vector< CohortIndividualSource::ColumnType > FromFileCohortIndividualSource::read_column_types( std::istream& stream, std::vector< std::string > const& column_names ) const {
		std::vector< ColumnType > result ;
		std::string line ;
		std::getline( stream, line ) ;
		std::vector< std::string > type_strings = string_utils::split_and_strip_discarding_empty_entries( line ) ;
		if( type_strings.size() != column_names.size() ) {
			throw MalformedInputError( m_filename, 1 ) ;
		}
		for( std::size_t i = 0; i < type_strings.size(); ++i ) {
			if( type_strings[i].size() != 1 ) {
				throw MalformedInputError( m_filename, 1, i ) ;
			}
			switch( type_strings[i][0] ) {
				case '0':
					result.push_back( e_ID_COLUMN ) ;
					break ;
				case 'D':
					result.push_back( e_DISCRETE_COVARIATE ) ;
					break ;
				case 'C':
					result.push_back( e_CONTINUOUS_COVARIATE ) ;
					break ;
				case 'B':
					result.push_back( e_BINARY_PHENOTYPE ) ;
					break ;
				case 'P':
					result.push_back( e_CONTINUOUS_PHENOTYPE ) ;
					break ;
				default:
					throw MalformedInputError( m_filename, 1, i ) ;
					break ;
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

	std::vector< std::vector< CohortIndividualSource::Entry > > FromFileCohortIndividualSource::read_entries( std::istream& stream, std::vector< ColumnType > const& column_types ) const {
		std::vector< std::vector< Entry > > result ;
		std::string line ;
		while( std::getline( stream, line ) ) {
			std::vector< std::string > string_entries = string_utils::split_and_strip_discarding_empty_entries( line ) ;
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
