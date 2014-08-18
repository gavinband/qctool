
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <set>
#include <fstream>
#include <boost/variant/variant.hpp>
#include <boost/bind.hpp>
#if HAVE_BOOST_SPIRIT
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/std_pair.hpp>
#endif
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
        std::size_t column_i = 0 ;
        if( genfile::string_utils::to_upper( column_name ) != "ID_1" ) {
            column_i = find_column_name( column_name ) ;
        }
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
		std::map< std::string, std::size_t >::const_iterator where = m_column_indices.find( string_utils::to_lower( column_name ) ) ;
		return ( where != m_column_indices.end() ) ;
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

	std::vector< std::size_t > FromFileCohortIndividualSource::find_samples_by_value( std::string const& column_name, Entry const& entry ) const {
        std::string const upperName = genfile::string_utils::to_upper( column_name ) ;
		if( upperName == "ID_1" || upperName == m_column_names[0] ) {
            // optimisation: fast lookup of principal ID field.
			std::map< Entry, std::size_t >::const_iterator where = m_sample_indices.find( entry ) ;
			if( where == m_sample_indices.end() ) {
				throw BadArgumentError( "FromFileCohortIndividualSource::find_entries()", "entry=\"" + string_utils::to_string( entry ) + "\"", "The entry was not found in column \"" + column_name + "\"" ) ;
			}
			return std::vector< std::size_t >( 1, where->second ) ;
		} else {
			return CohortIndividualSource::find_samples_by_value( column_name, entry ) ;
		}
	}

	std::size_t FromFileCohortIndividualSource::find_column_name( std::string const& column_name ) const {
		std::map< std::string, std::size_t >::const_iterator where = m_column_indices.find( string_utils::to_lower( column_name ) ) ;
		if( where == m_column_indices.end() ) {
			throw BadArgumentError( "FromFileCohortIndividualSource::find_column_name()", "column_name=\"" + column_name + "\"", "Column \"" + column_name + "\" not found" ) ;
		}
		return where->second ;
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
		// read any header comments
		m_comments = read_comments( stream ) ;
		m_number_of_metadata_lines = std::count( m_comments.begin(), m_comments.end(), '\n' ) ;
		
		typedef std::map< std::string, CohortIndividualSource::ColumnType > ColumnTypeMap ;
		boost::optional< ColumnTypeMap > column_types = read_column_types_from_comments( m_comments ) ;
		
		m_column_names = read_column_names( stream ) ;
		for( std::size_t i = 0; i < m_column_names.size(); ++i ) {
			if( !m_column_indices.insert( std::make_pair( string_utils::to_lower( m_column_names[i] ), i ) ).second ) {
				throw MalformedInputError( m_filename, "Column \"" + string_utils::to_lower( m_column_names[i] ) + "\" appears more than once", m_number_of_metadata_lines + 0 ) ;
			}
		}

		if( !column_types ) {
			column_types = read_column_type_line( stream, m_column_names ) ;
		}

		if( !column_types ) {
			throw MalformedInputError( m_filename, "Expected column type line after column names.", m_number_of_metadata_lines + 1 ) ;
        }

		m_column_types.resize( m_column_names.size() ) ;
		for( std::size_t i = 0; i < m_column_names.size(); ++i ) {
			ColumnTypeMap::const_iterator where = (*column_types).find( m_column_names[i] ) ;
			if( where == (*column_types).end() ) {
				throw MalformedInputError( m_filename, m_number_of_metadata_lines + 1 ) ;
			}
			m_column_types[i] = where->second ;
			
            if( genfile::string_utils::to_upper( m_column_names[i] ) == "ID_1" && i != 0 ) {
				throw MalformedInputError( m_filename, "The name 'ID_1', if used, must be the first column in the sample file.", 0 ) ;
            }
            if( genfile::string_utils::to_upper( m_column_names[i] ) == "ID_2" && i != 1 ) {
				throw MalformedInputError( m_filename, "The name 'ID_2', if used, must be the second column in the sample file.", 0 ) ;
            }
            if( genfile::string_utils::to_upper( m_column_names[i] ) == "missing" && i != 2 ) {
				throw MalformedInputError( m_filename, "The name 'missing', if used, must be the third column in the sample file.", 0 ) ;
            }
			if( i < 3 ) {
				if( m_column_types[i] != e_ID_COLUMN ) {
					throw MalformedInputError( m_filename, "column " + genfile::string_utils::to_string( i+1 ) + "(\"" + m_column_names[i] + "\") should have type \"0\"", 0 ) ;
				}
			} else {
				if( m_column_types[i] == e_ID_COLUMN ) {
					throw MalformedInputError( m_filename, "Only the first (ID_1), second (ID_2), and third (missing) columns can have type \"0\"", 0 ) ;
				}
			}
            // update the missingness column to the right type.
			if( genfile::string_utils::to_upper( m_column_names[i] ) == "MISSING" ) {
				m_column_types[i] = e_MISSINGNESS_COLUMN ;
			}
		}

		assert( m_column_names.size() == m_column_types.size() ) ;
		m_entries = read_entries( stream, m_column_types ) ;
		for( std::size_t i = 0; i < m_entries.size(); ++i ) {
			m_sample_indices[ m_entries[i].front() ] = i ;
		}
		assert( stream.eof() ) ;
	}
	
	std::string FromFileCohortIndividualSource::read_comments( std::istream& stream ) const {
		std::string result ;
		while( stream.peek() == '#' ) {
			std::string line ;
			std::getline( stream, line ) ;
			assert( line.size() > 0 ) ;
			result.append( line.substr( 1, line.size() ) + "\n" ) ;
		}
		return result ;
	}
	
	namespace {
		boost::optional< CohortIndividualSource::ColumnType > get_column_type( std::string const& type_string ) {
			boost::optional< CohortIndividualSource::ColumnType > result ;
			if( type_string.size() == 1 ) {
				switch( type_string[0] ) {
					case '0':
						result = CohortIndividualSource::e_ID_COLUMN ;
						break ;
					case 'D':
						result = CohortIndividualSource::e_DISCRETE_COVARIATE ;
						break ;
					case 'C':
						result = CohortIndividualSource::e_CONTINUOUS_COVARIATE ;
						break ;
					case 'B':
						result = CohortIndividualSource::e_BINARY_PHENOTYPE ;
						break ;
					case 'P':
						result = CohortIndividualSource::e_CONTINUOUS_PHENOTYPE ;
						break ;
					default:
						break ;
				}
			}
			return result ;
		}

		void insert_name_and_type(
			std::map< std::string, CohortIndividualSource::ColumnType >* result,
			std::vector< char > const& the_name,
			std::vector< char > const& the_type
		) {
			assert( result ) ;
			std::string name( the_name.begin(), the_name.end() ) ;
			std::string type_string( the_type.begin(), the_type.end() ) ;
			std::cerr << "name:" << name << ", type:" << type_string << ".\n" ;
			
			boost::optional< CohortIndividualSource::ColumnType > type = get_column_type( type_string ) ;
			if( !type ) {
				throw MalformedInputError( "(unknown)", 0 ) ;
			}
			(*result)[ name ] = (*type) ;
		}
	}

	boost::optional< std::map< std::string, CohortIndividualSource::ColumnType > > FromFileCohortIndividualSource::read_column_types_from_comments(
		std::string const& comments
	) const {
		boost::optional< std::map< std::string, CohortIndividualSource::ColumnType > > result ;
#if !HAVE_BOOST_SPIRIT
		return result ;
#else
		using namespace boost::spirit::qi ;
		using namespace boost::phoenix ;
		using boost::spirit::ascii::char_ ;
		using boost::spirit::ascii::graph ;
		using boost::spirit::ascii::blank ;

		typedef std::map< std::string, CohortIndividualSource::ColumnType > Intermediate ;

		std::size_t pos = comments.find( "metadata" ) ;
		if( pos != std::string::npos ) {
			Intermediate intermediate ;
			std::string::const_iterator begin = comments.begin() + pos ;
			bool parsed = boost::spirit::qi::phrase_parse(
				begin, comments.end(),
					lit( "metadata:" )
					>> '{'
					>> (
						lit( "\"" )
						>> ( +(graph - '\"') ) >> '\"'
						>> ':'
						>> '{'
						>> "\"type\""
						>> ':'
						>> '\"'
						>> ( +(graph - '\"'))
						>> '\"'
						>> '}'
					)
					[
						boost::phoenix::bind(
							&insert_name_and_type,
							&intermediate,
							boost::spirit::_1,
							boost::spirit::_2
						) 
					]
					% ','
				>> '}',
				space
			) ;
	
			if( parsed ) {
				result = intermediate ;
			} else {
				throw MalformedInputError( m_filename, "\"metadata\" line found in header, but column specification was malformed", 0 ) ;
			}
		}
		return result ;
#endif
	}
	
	std::vector< std::string > FromFileCohortIndividualSource::read_column_names( std::istream& stream ) const {
		std::vector< std::string > result ;
		std::string line ;
		std::getline( stream, line ) ;
		result = string_utils::split_and_strip_discarding_empty_entries( line ) ;
		if( result.size() < 1 ) {
			throw MalformedInputError( m_filename, 0 + m_number_of_metadata_lines ) ;
		}
		return result ;
	}

	std::map< std::string, CohortIndividualSource::ColumnType > FromFileCohortIndividualSource::read_column_type_line( std::istream& stream, std::vector< std::string > const& column_names ) const {
		std::map< std::string, CohortIndividualSource::ColumnType > result ;
		std::string line ;
		std::getline( stream, line ) ;
		std::vector< std::string > type_strings = string_utils::split_and_strip_discarding_empty_entries( line ) ;
		if( type_strings.size() != column_names.size() ) {
			throw MalformedInputError( m_filename, 1 + m_number_of_metadata_lines ) ;
		}
		for( std::size_t i = 0; i < type_strings.size(); ++i ) {
			boost::optional< CohortIndividualSource::ColumnType > type = get_column_type( type_strings[i] ) ;
			if( !type ) {
				throw MalformedInputError( m_filename, 1 + m_number_of_metadata_lines, i ) ;
			} else {
				result[ m_column_names[i] ] = *type ;
			}
		}
		assert( result.size() == column_names.size() ) ;

		return result ;
	}

	std::vector< std::vector< CohortIndividualSource::Entry > > FromFileCohortIndividualSource::read_entries( std::istream& stream, std::vector< ColumnType > const& column_types ) const {
		std::vector< std::vector< Entry > > result ;
		std::string line ;
		while( std::getline( stream, line ) ) {
			std::vector< std::string > string_entries = string_utils::split_and_strip_discarding_empty_entries( line ) ;
			if( string_entries.size() != column_types.size() ) {
				throw MalformedInputError( m_filename, 2 + m_number_of_metadata_lines + result.size() ) ;
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
				throw MalformedInputError( m_filename, line_number + m_number_of_metadata_lines, i ) ;
			}
			catch( UnexpectedMissingValueError const& e ) {
				throw UnexpectedMissingValueError( m_filename, line_number + m_number_of_metadata_lines, i ) ;
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
