#include <string>
#include <map>

#include <memory>
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/snp_data_utils.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"
#include "genfile/VCFFormatMetaDataParser.hpp"

using boost::tuples::tie ;

namespace genfile {
	VCFFormatMetaDataParser::VCFFormatMetaDataParser( std::string const& spec, std::istream& stream ):
		m_spec( spec ),
		m_metadata( read_metadata( stream ))
	{}
	
	VCFFormatMetaDataParser::Metadata VCFFormatMetaDataParser::read_metadata( std::istream& in ) const {
		std::istream::iostate old_exceptions = in.exceptions() ;
		in.exceptions( std::ios::eofbit | std::ios::failbit | std::ios::badbit ) ;
		Metadata result ;
		
		for( std::size_t line_number = 0; in; ++line_number ) {
			try {
				char a_char = static_cast< char >( in.get() ) ;
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
					result.insert( parse_meta_line( line_number, line )) ;
				} else {
					in.putback( a_char ) ;
					break ;
				}
			}
			catch( std::exception const& ) {
				throw MalformedInputError( m_spec, line_number ) ;
			}
		}
		
		if( result.empty() ) {
			throw MalformedInputError( m_spec, 0 ) ;
		}
		
		in.exceptions( old_exceptions ) ;
		return result ;
	}
	
	std::pair< std::string, std::map< std::string, std::string > > VCFFormatMetaDataParser::parse_meta_line(
		std::size_t line_number,
		std::string const& line
	) const {
		std::string key, value ;
		tie( key, value ) = parse_key_value_pair( line ) ;
		if( !validate_meta_key( line_number, key )) {
			throw MalformedInputError( m_spec, line_number ) ;
		} ;
		
		std::map< std::string, std::string > result ;
		
		if( key == "fileformat" ) {
			if( value != "VCFv4.1" ) {
				throw FormatUnsupportedError( m_spec, value ) ;
			}
		}

		if( key == "INFO" || key == "FILTER" || key == "FORMAT") {
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

	std::pair< std::string, std::string > VCFFormatMetaDataParser::parse_key_value_pair( std::string const& line ) const {
		std::size_t pos = line.find( '=' ) ;
		if( pos == std::string::npos ) {
			throw BadArgumentError( "VCFFormatMetaDataParser::parse_key_value_pair()", "line = \"" + line + "\"" ) ;
		}
		return std::make_pair( line.substr( 0, pos ), line.substr( pos + 1, line.size() )) ;
	}

	bool VCFFormatMetaDataParser::validate_meta_key(
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
	
	std::map< std::string, std::string > VCFFormatMetaDataParser::parse_meta_value(
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
				result.insert( parse_key_value_pair( elts[i] ) ) ;
			}
		}
		catch( BadArgumentError const& e ) {
			throw MalformedInputError( m_spec, line_number ) ;
		}
		return result ;
	}
	
	bool VCFFormatMetaDataParser::validate_meta_value(
		std::string const& key,
		std::map< std::string,
		std::string > const& meta_value
	) const {
		std::map< std::string, std::string >::const_iterator where ;
		// validate ID
		if( ( where = meta_value.find( "ID" ) ) == meta_value.end() ) { return false ; }
		// validate Description
		if( ( where = meta_value.find( "Description" ) ) == meta_value.end() ) { return false ; }
		if( where->second.size() < 2 || where->second[0] != '"' || where->second[ where->second.size() - 1 ] != '"' ) {
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
			if( ( where = meta_value.find( "Number" ) ) == meta_value.end() ) { return false ; }
			if( where->second != "A" && where->second != "G" && where->second != "." ) {
				// Must be an integer.
				try {
					string_utils::to_repr< unsigned int >( where->second ) ;
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
					|| ( key == "INFO" && where->second == "Flag" )
				)
			) {
				return false ;
			}
		}
		return true ;
	}
}
