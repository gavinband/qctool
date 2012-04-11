
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/vcf/TrivialMetadataParser.hpp"
#include "genfile/FileUtils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	namespace vcf {
		TrivialMetadataParser::TrivialMetadataParser( std::string const& filename ):
			m_spec( filename )
		{
			std::auto_ptr< std::istream > file = open_text_file_for_input( filename ) ;
			setup( *file ) ;
		}
		
		TrivialMetadataParser::TrivialMetadataParser( std::string const& spec, std::istream& stream ):
			m_spec( spec )
		{
			setup( stream ) ;
		}

		void TrivialMetadataParser::setup( std::istream& stream ) {
			std::istream::iostate old_exceptions = stream.exceptions() ;
			for( m_number_of_lines = 0; stream; ++m_number_of_lines ) {
				if( !read_metadata_line( stream, m_number_of_lines ) ) {
					break ;
				}
			}
			stream.exceptions( old_exceptions ) ;
		}

		bool TrivialMetadataParser::read_metadata_line( std::istream& in, std::size_t line_number ) const {
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
	}
}
