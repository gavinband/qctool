#include "genfile/StringView.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <cassert>

namespace genfile {
	namespace stringview {
		StringView::StringView( std::string const& a_string ):
			m_string( &a_string ),
			m_start( 0 ),
			m_end( a_string.size() )
		{}
		
		StringView::StringView( std::string const& a_string, std::size_t start, std::size_t end ):
			m_string( &a_string ),
			m_start( start ),
			m_end( end )
		{
			assert( start <= end ) ;
			assert( end <= m_string->size() ) ;
		}

		StringView::StringView( StringView const& other, std::size_t start, std::size_t end ):
			m_string( other.m_string ),
			m_start( start + other.m_start ),
			m_end( end + other.m_start)
		{
			assert( start <= end ) ;
			assert( end <= m_string->size() ) ;
		}

		StringView::StringView( StringView const& other ):
			m_string( other.m_string ),
			m_start( other.m_start ),
			m_end( other.m_end )
		{}
		
		StringView& StringView::operator=( StringView const& other ) {
			m_string = other.m_string ;
			m_start = other.m_start ;
			m_end = other.m_end ;
			return *this ;
		}
		
		std::size_t StringView::find( char c, std::size_t pos ) const {
			for( std::size_t i = m_start + pos; i < m_end; ++i ) {
				if( (*m_string)[i] == c ) {
					return i - m_start ;
				}
			}
			return std::string::npos ;
		}

		std::size_t StringView::find_first_of( std::string const& chars, std::size_t pos ) const {
			unsigned char array[ 256 ] ;
			std::memset( array, 0, 256 ) ;
			for( char const* i = &chars[0]; i < &chars[0] + chars.size(); ++i ) {
				array[ static_cast< int >( *i ) ] = 1 ;
			}
			for( std::size_t i = m_start + pos; i < m_end; ++i ) {
				if( array[ static_cast< int >( (*m_string)[i] ) ] ) {
					return i - m_start ;
				}
			}
			return std::string::npos ;
		}
		
		bool StringView::operator==( std::string const& other ) const {
			if( size() != other.size() ) {
				return false ;
			}
			for( std::size_t i = m_start; i < m_end; ++i ) {
				if( (*m_string)[i] != other[ i - m_start ] ) {
					return false ;
				}
			}
			return true ;
		}
		
		std::vector< StringView > split( StringView const& string_to_split, std::string const& split_chars ) {
			std::vector< StringView > result ;
			std::size_t last_pos = 0, pos = 0 ;
			for( ; pos != string_to_split.size(); last_pos = pos + 1 ) {
				pos = string_to_split.find_first_of( split_chars, last_pos ) ;
				if( pos == std::string::npos ) {
					pos = string_to_split.size() ;
				}
				result.push_back( StringView( string_to_split, last_pos, pos ) ) ;
			}

			return result ;	
		}
		
		
		
	}
}
