
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <map>

#include <memory>
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/snp_data_utils.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"
#include "genfile/vcf/StrictMetadataParser.hpp"

using boost::tuples::tie ;

namespace genfile {
	namespace vcf {
		StrictMetadataParser::StrictMetadataParser( std::string const& filename ):
			m_spec( filename )
		{
			std::auto_ptr< std::istream > stream = open_text_file_for_input( filename ) ;
			setup( filename, *stream ) ;
		}

		StrictMetadataParser::StrictMetadataParser( std::string const& spec, std::istream& stream ):
			m_spec( spec )
		{
			setup( spec, stream ) ;
		}
		
		void StrictMetadataParser::setup( std::string const& spec, std::istream& stream ) {
			assert( stream ) ;
			m_version = read_version( stream ) ;
			m_metadata = read_metadata( stream, m_version ) ;
		}
	
		std::string StrictMetadataParser::read_version( std::istream& in ) const {
			std::istream::iostate old_exceptions = in.exceptions() ;
			in.exceptions( std::ios::eofbit | std::ios::failbit | std::ios::badbit ) ;
			Metadata fileformat_spec ;
			if( !read_metadata_line( in, 0, &fileformat_spec )) {
				throw MalformedInputError( m_spec, 0 ) ;
			}
			assert( fileformat_spec.size() == 1 ) ;
			Metadata::const_iterator i = fileformat_spec.find( "fileformat" ) ;
			if( i == fileformat_spec.end() ) {
				throw MalformedInputError( m_spec, 0 ) ;
			}
			std::map< std::string, std::string >::const_iterator j = i->second.find( "version" ) ;
			if( j == i->second.end() ) {
				throw MalformedInputError( m_spec, 0 ) ;
			}
			std::string const result = j->second ;
			assert( result == "4.0" || result == "4.1" || result == "4.2" ) ;
			double const version = string_utils::to_repr< double >( result ) ;
			if( version < 4.0 ) {
				throw MalformedInputError( m_spec, 0 ) ;
			}
			in.exceptions( old_exceptions ) ;
			return result ;
		}

		StrictMetadataParser::Metadata StrictMetadataParser::read_metadata( std::istream& in, std::string const& version ) const {
			std::istream::iostate old_exceptions = in.exceptions() ;
			Metadata result ;
			{
				// put version information into the result.
				std::map< std::string, std::string > A ;
				A[ "version" ] = version ;
				result.insert( std::make_pair( "fileformat", A )) ;
			}
			for( std::size_t line_number = 1; in; ++line_number ) {
				if( !read_metadata_line( in, line_number, &result ) ) {
					break ;
				}
			}
			if( result.empty() ) {
				throw MalformedInputError( m_spec, 0 ) ;
			}
		
			in.exceptions( old_exceptions ) ;
			return result ;
		}

		bool StrictMetadataParser::read_metadata_line( std::istream& in, std::size_t line_number, Metadata* result ) const {
			try {
				in.exceptions( std::ios::badbit ) ;
				char a_char = static_cast< char >( in.get() ) ;
				if( !in ) {
					return false ;
				}
				in.exceptions( std::ios::eofbit | std::ios::failbit | std::ios::badbit ) ;
				
				char next_char = static_cast< char >( in.peek() ) ;
				if( a_char == '\t' || next_char == '\t' ) {
					throw MalformedInputError( m_spec, line_number ) ;
				}
				if( a_char == '#' && next_char == '#' ) {
					in.get() ;
					std::string line ;
					std::getline( in, line ) ;
					if( line.find( "\t" ) != std::string::npos ) {
						throw MalformedInputError( m_spec, line_number ) ;
					}
					try {
						result->insert( parse_meta_line( line_number, line )) ;
					}
					catch( MalformedInputError const& error ) {
						std::map< std::string, std::string > data ;
						data[ "what" ] = line ;
						result->insert( std::make_pair( "parse_error", data ) ) ;
						// throw ;
					}
				} else {
					in.putback( a_char ) ;
					return false ;
				}
			}
			catch( std::ios_base::failure const& ) {
				throw MalformedInputError( m_spec, line_number ) ;
			}
			return true ;
		}
	
		std::pair< std::string, std::map< std::string, std::string > > StrictMetadataParser::parse_meta_line(
			std::size_t line_number,
			std::string const& line
		) const {
			std::string key, value ;
			tie( key, value ) = parse_key_value_pair( line_number, line ) ;
			if( !validate_meta_key( line_number, key )) {
				throw MalformedInputError( m_spec, line_number ) ;
			} else {
				std::map< std::string, std::string > result ;
			
				if( key == "fileformat" ) {
					if( value != "VCFv4.2" && value != "VCFv4.1" && value != "VCFv4.0" ) {
						throw FormatUnsupportedError( m_spec, value ) ;
					}
					result[ "version" ] = value.substr( 4, 3 ) ;
				}
				else if( key == "INFO" || key == "FILTER" || key == "FORMAT") {
					result = parse_meta_value( line_number, key, value ) ;
					if( !validate_meta_value( key, result )) {
						throw MalformedInputError( m_spec, line_number ) ;
					}
				}
				else {
					// If we get here, unknown tag type.
					// Return a map with key "" and value the whole string.
					result[ "" ] = value ;
				}
				return std::make_pair( key, result ) ;
			}
		}

		std::pair< std::string, std::string > StrictMetadataParser::parse_key_value_pair( std::size_t const line_number, std::string const& line ) const {
			std::size_t pos = line.find( '=' ) ;
			if( pos == std::string::npos ) {
				// Take "" as key, whole line as value.
				return std::make_pair( "", line ) ;
			}
			std::pair< std::string, std::string > result = std::make_pair( line.substr( 0, pos ), line.substr( pos + 1, line.size() )) ;
			if( result.first == "Description" ) {
				if( result.second.size() < 2 || result.second[0] != '"' || result.second[result.second.size() - 1] != '"' ) {
					std::cerr << "tt\n" ;
					throw MalformedInputError( m_spec, line_number ) ;
				}
				result.second = result.second.substr( 1, result.second.size() - 2 ) ;
			}
			return result ;
		}

		bool StrictMetadataParser::validate_meta_key(
			std::size_t line_number,
			std::string const& key
		) const {
			if( key.empty() ) {
				return false ;
			}
			else if( line_number == 0 && key != "fileformat" ) {
				return false ;
			}
			return true ;
		}
	
		std::map< std::string, std::string > StrictMetadataParser::parse_meta_value(
			std::size_t line_number,
			std::string const& key,
			std::string value
		) const {
			std::map< std::string, std::string > result ;
			if( value.size() < 2 || value[0] != '<' || value[ value.size() - 1 ] != '>' ) {
				throw MalformedInputError( m_spec, line_number ) ;
			}
			try {
				// strip angle brackets.
				std::vector< std::string > elts = string_utils::split_respecting_delimited_regions(
					value.substr( 1, value.size() - 2 ),
					",",
					"\"\""
				) ;
				for( std::size_t i = 0; i < elts.size(); ++i ) {
					result.insert( parse_key_value_pair( line_number, elts[i] ) ) ;
				}
			}
			catch( BadArgumentError const& e ) {
				throw MalformedInputError( m_spec, line_number ) ;
			}
			return result ;
		}
	
		bool StrictMetadataParser::validate_meta_value(
			std::string const& key,
			std::map< std::string,
			std::string > const& meta_value
		) const {
			std::map< std::string, std::string >::const_iterator where ;
			// validate ID
			if( ( where = meta_value.find( "ID" ) ) == meta_value.end() ) { return false ; }
			// validate Description
			if( ( where = meta_value.find( "Description" ) ) == meta_value.end() ) { return false ; }
			if( where->second.size() > 1 && ( where->second[0] == '"' || where->second[ where->second.size() - 1 ] == '"' )) {
				return false ;
			}
			if( key == "FILTER" ) {
				if( meta_value.size() != 2 ) {
					return false ;
				}
			}
			else if( key == "INFO" || key == "FORMAT" ) {
				if( meta_value.size() != 4 ) {
					return false ;
				}
				// validate Number
				if( ( where = meta_value.find( "Number" ) ) == meta_value.end() ) {
					return false ;
				}
				if(
					where->second == "."
					|| ( ( m_version == "4.1" || m_version == "4.2" ) && ( where->second == "A" || where->second == "G" ))
					|| ( ( m_version == "4.2" ) && ( where->second == "R" ))
				) {
					// Ok, a non-numerical type.
				}
				else {
					// Must be an integer.
					try {
						int const value = string_utils::to_repr< int >( where->second ) ;
						if( value < 0 ) { return false ; }
					} catch( string_utils::StringConversionError const& ) {
						return false ;
					}
				}
				// validate Type
				if( ( where = meta_value.find( "Type" ) ) == meta_value.end() ) { return false ; }
				if( !(
						where->second == "Integer"
						|| where->second == "Float"
						|| where->second == "Character"
						|| where->second == "String"
						|| where->second == "Genotype"
						|| ( key == "INFO" && where->second == "Flag" )
					)
				) {
					return false ;
				}
			}
			return true ;
		}
	}
}
