#include <cstring>
#include <limits>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <boost/bind.hpp>
#include "genfile/string_utils/slice.hpp"


namespace genfile {
	namespace string_utils {
		slice::slice( std::string const& a_string ):
			m_string( &a_string ),
			m_start( 0 ),
			m_end( a_string.size() )
		{}
		
		slice::slice( std::string const& a_string, std::size_t start, std::size_t end ):
			m_string( &a_string ),
			m_start( start ),
			m_end( end )
		{
			assert( start <= end ) ;
			assert( end <= m_string->size() ) ;
		}

		slice::slice( slice const& other, std::size_t start, std::size_t end ):
			m_string( other.m_string ),
			m_start( start + other.m_start ),
			m_end( end + other.m_start )
		{
			assert( start <= end ) ;
			assert( end <= m_string->size() ) ;
		}

		slice::slice( slice const& other ):
			m_string( other.m_string ),
			m_start( other.m_start ),
			m_end( other.m_end )
		{}
		
		/*
		slice& slice::operator=( slice const& other ) {
			m_string = other.m_string ;
			m_start = other.m_start ;
			m_end = other.m_end ;
			return *this ;
		}
		*/
		
		std::size_t slice::find( char c, std::size_t pos ) const {
			for( std::size_t i = m_start + pos; i < m_end; ++i ) {
				if( (*m_string)[i] == c ) {
					return i - m_start ;
				}
			}
			return std::string::npos ;
		}

		namespace impl {
			void make_membership_array( std::string const& chars, char* array, std::size_t const n ) {
				assert( n >= std::size_t( std::numeric_limits< unsigned char >::max() + 1 )) ;
				std::memset( array, 0, n ) ;
				for( char const* i = &chars[0]; i < (&chars[0] + chars.size()); ++i ) {
					array[ static_cast< int >( static_cast< unsigned char const >( *i ) ) ] = 1 ;
				}
			}
		}

		std::size_t slice::find_first_of( std::string const& chars, std::size_t pos ) const {
			assert( std::numeric_limits< unsigned char >::max() + 1 == 256 ) ;
			char array[ 256 ] ;
			impl::make_membership_array( chars, array, 256 ) ;
			return find_first_of( array, pos ) ;
		}
		
		std::size_t slice::find_first_of( char* membership_array, std::size_t pos ) const {
			for( std::size_t i = m_start + pos; i < m_end; ++i ) {
				if( membership_array[ static_cast< int >( (*m_string)[i] ) ] ) {
					return i - m_start ;
				}
			}
			return std::string::npos ;
		}

		bool slice::operator==( std::string const& other ) const {
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

		bool slice::operator!=( std::string const& other ) const {
			return !( *this == other ) ;
		}

		bool slice::operator==( slice const& other ) const {
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

		std::vector< slice > slice::split( std::string const& split_chars ) const {
			std::vector< slice > result ;
			result.reserve( 1000 ) ;
			assert( std::numeric_limits< unsigned char >::max() + 1 == 256 ) ;
			char array[ 256 ] ;
			impl::make_membership_array( split_chars, array, 256 ) ;
			std::size_t last_pos = 0, pos = 0 ;
			do {
				pos = find_first_of( array, last_pos ) ;
				if( pos == std::string::npos ) {
					pos = size() ;
				}
				result.push_back( slice( *this, last_pos, pos ) ) ;
				last_pos = pos + 1 ;
			}
			while( pos != size() ) ;
			return result ;	
		}

		void slice::split( std::string const& split_chars, std::vector< slice >* result ) const {
			assert( std::numeric_limits< unsigned char >::max() + 1 == 256 ) ;
			char array[ 256 ] ;
			impl::make_membership_array( split_chars, array, 256 ) ;
			std::size_t last_pos = 0, pos = 0 ;
			do {
				pos = find_first_of( array, last_pos ) ;
				if( pos == std::string::npos ) {
					pos = size() ;
				}
				result->push_back( slice( *this, last_pos, pos ) ) ;
				last_pos = pos + 1 ;
			}
			while( pos != size() ) ;
		}

		
		
		
	}
}
