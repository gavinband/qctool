
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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
		
		slice::slice( char const* c_string ):
			m_data( c_string ),
			m_end_of_data( c_string + std::strlen( c_string ) ),
			m_start( 0 )
		{
			assert( c_string ) ;
			m_end = (m_end_of_data - m_data) ;
		}
		
		slice::slice( char const* c_string, std::size_t start, std::size_t end ):
			m_data( c_string ),
			m_end_of_data( c_string + std::strlen( c_string ) ),
			m_start( start ),
			m_end( end )
		{
			assert( c_string ) ;
			assert( m_start <= m_end ) ;
			assert( m_data + m_end <= m_end_of_data ) ;
		}
		
		slice::slice( char const* data, char const* end ):
			m_data( data ),
			m_end_of_data( end ),
			m_start( 0 ),
			m_end( end - m_data )
		{
			assert( data ) ;
			assert( end ) ;
			assert( data <= end ) ;
		}
			
		// A view of a subset of a character buffer
		slice::slice( char const* data, char const* end_of_data, std::size_t start, std::size_t end ):
			m_data( data ),
			m_end_of_data( end_of_data ),
			m_start( start ),
			m_end( end )
		{
			assert( data ) ;
			assert( end_of_data ) ;
			assert( data <= end_of_data ) ;
			assert( start <= end ) ;
			assert( data + end <= end_of_data ) ;
		}

		slice::slice( std::string const& a_string ):
			m_data( a_string.data() ),
			m_end_of_data( a_string.data() + a_string.size() ),
			m_start( 0 ),
			m_end( a_string.size() )
		{}
		
		slice::slice( std::string const& a_string, std::size_t start, std::size_t end ):
			m_data( a_string.data() ),
			m_end_of_data( a_string.data() + a_string.size() ),
			m_start( start ),
			m_end( end )
		{
			assert( start <= end ) ;
			assert( m_data + m_end <= m_end_of_data ) ;
		}

		slice::slice( slice const& other, std::size_t start, std::size_t end ):
			m_data( other.m_data ),
			m_end_of_data( other.m_end_of_data ),
			m_start( start + other.m_start ),
			m_end( end + other.m_start )
		{
			assert( start <= end ) ;
			assert( m_data + m_end <= m_end_of_data ) ;
		}

		slice::slice( slice const& other ):
			m_data( other.m_data ),
			m_end_of_data( other.m_end_of_data ),
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
				if( *(m_data + i) == c ) {
					return i - m_start ;
				}
			}
			return std::string::npos ;
		}

		namespace impl {
			void make_membership_array( std::string const& chars, char* array, std::size_t const n, char const out_value = 0, char const in_value = 1 ) {
				assert( n >= std::size_t( std::numeric_limits< unsigned char >::max() + 1 )) ;
				std::memset( array, out_value, n ) ;
				for( char const* i = &chars[0]; i < (&chars[0] + chars.size()); ++i ) {
					array[ static_cast< int >( static_cast< unsigned char const >( *i ) ) ] = in_value ;
				}
			}
		}

		std::size_t slice::find_first_of( std::string const& chars, std::size_t pos ) const {
			assert( std::numeric_limits< unsigned char >::max() + 1 == 256 ) ;
			char array[ 256 ] ;
			impl::make_membership_array( chars, array, 256, 0, 1 ) ;
			return find_first_of( array, pos ) ;
		}

		std::size_t slice::find_first_not_of( std::string const& chars, std::size_t pos ) const {
			assert( std::numeric_limits< unsigned char >::max() + 1 == 256 ) ;
			char array[ 256 ] ;
			impl::make_membership_array( chars, array, 256, 1, 0 ) ;
			return find_first_of( array, pos ) ;
		}
		
		std::size_t slice::find_first_of( char* membership_array, std::size_t pos ) const {
			for( std::size_t i = m_start + pos; i < m_end; ++i ) {
				if( membership_array[ static_cast< int >( *(m_data + i) ) ] ) {
					return i - m_start ;
				}
			}
			return std::string::npos ;
		}

		std::size_t slice::find_last_of( std::string const& chars, std::size_t pos ) const {
			assert( std::numeric_limits< unsigned char >::max() + 1 == 256 ) ;
			char array[ 256 ] ;
			impl::make_membership_array( chars, array, 256, 0, 1 ) ;
			return find_last_of( array, pos ) ;
		}

		std::size_t slice::find_last_not_of( std::string const& chars, std::size_t pos ) const {
			assert( std::numeric_limits< unsigned char >::max() + 1 == 256 ) ;
			char array[ 256 ] ;
			impl::make_membership_array( chars, array, 256, 1, 0 ) ;
			return find_last_of( array, pos ) ;
		}
		
		std::size_t slice::find_last_of( char* membership_array, std::size_t pos ) const {
			if( pos == std::string::npos ) {
				pos = m_end - m_start - 1 ;
			}
			for( std::size_t i = m_start + pos + 1; i > m_start; --i ) {
				if( membership_array[ static_cast< int >( *(m_data+i-1) ) ] ) {
					return i - 1 - m_start ;
				}
			}
			return std::string::npos ;
		}

		bool slice::operator==( std::string const& other ) const {
			if( size() != other.size() ) {
				return false ;
			}
			for( std::size_t i = m_start; i < m_end; ++i ) {
				if( *(m_data+i) != other[ i - m_start ] ) {
					return false ;
				}
			}
			return true ;
		}

		bool slice::operator!=( std::string const& other ) const {
			return !( *this == other ) ;
		}

		bool slice::operator==( char const* other ) const {
			int const other_length = std::strlen( other ) ;
			if( size() != other_length ) {
				return false ;
			}
			for( std::size_t i = m_start; i < m_end; ++i ) {
				if( *(m_data+i) != *(other+i-m_start) ) {
					return false ;
				}
			}
			return true ;
		}

		bool slice::operator!=( char const* other ) const {
			return !( *this == other ) ;
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

		slice slice::substr( std::size_t start, std::size_t end ) const {
			return slice( *this, start, end ) ;
		}

		slice slice::strip( std::string const& chars ) const {
			if( chars.empty() ) {
				return *this ;
			}

			std::size_t lpos = find_first_not_of( chars ) ;
			if( lpos == std::string::npos ) {
				return slice( m_data, m_end_of_data, m_start, m_start ) ;
			}

			std::size_t rpos = find_last_not_of( chars ) ;
			return slice( m_data, m_end_of_data, m_start + lpos, m_start + rpos + 1 ) ;
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
		
		void slice::split( std::string const& split_chars, boost::function< void( slice ) > callback ) const {
			assert( std::numeric_limits< unsigned char >::max() + 1 == 256 ) ;
			char array[ 256 ] ;
			impl::make_membership_array( split_chars, array, 256 ) ;
			std::size_t last_pos = 0, pos = 0 ;
			do {
				pos = find_first_of( array, last_pos ) ;
				if( pos == std::string::npos ) {
					pos = size() ;
				}
				callback( slice( *this, last_pos, pos ) ) ;
				last_pos = pos + 1 ;
			}
			while( pos != size() ) ;
		}

		slice::const_iterator slice::begin() const {
			return m_data + m_start ;
		}

		slice::const_iterator slice::end() const {
			return m_data + m_end ;
		}
		
		bool operator==( slice const& left, slice const& right ) {
			return (left.size() == right.size())
				&& (
					std::memcmp(
						left.m_data + left.m_start,
						right.m_data + right.m_start,
						left.size()
					) == 0
				) ;
		}

		bool operator!=( slice const& left, slice const& right ) {
			return (left.size() != right.size())
				|| (
					std::memcmp(
						left.m_data + left.m_start,
						right.m_data + right.m_start,
						left.size()
					) != 0
				) ;
		}

		bool operator<( slice const& left, slice const& right ) {
			int cmp = std::memcmp(
					left.m_data + left.m_start,
					right.m_data + right.m_start,
					std::min( left.size(), right.size() )
				) ;
			return (cmp < 0) || (cmp == 0 && left.size() < right.size()) ;
		}

		bool operator>( slice const& left, slice const& right ) {
			int cmp = std::memcmp(
					left.m_data + left.m_start,
					right.m_data + right.m_start,
					std::min( left.size(), right.size() )
				) ;
			return (cmp > 0) || (cmp == 0 && left.size() > right.size()) ;
		}

		bool operator<=( slice const& left, slice const& right ) {
			int cmp = std::memcmp(
					left.m_data + left.m_start,
					right.m_data + right.m_start,
					std::min( left.size(), right.size() )
				) ;
			return (cmp < 0) || (cmp == 0 && left.size() <= right.size()) ;
		}

		bool operator>=( slice const& left, slice const& right ) {
			int cmp = std::memcmp(
					left.m_data + left.m_start,
					right.m_data + right.m_start,
					std::min( left.size(), right.size() )
				) ;
			return (cmp > 0) || (cmp == 0 && left.size() >= right.size()) ;
		}
		
		std::ostream& operator<<( std::ostream& o, slice const& s ) {
			for( std::size_t i = s.m_start; i < s.m_end; ++i ) {
				o << (*(s.m_data+i)) ;
			}
			return o ;
		}
		
		std::string join( std::vector< slice > const& slices, std::string const& joiner ) {
			std::size_t size = 0 ;
			for( std::size_t i = 0; i < slices.size(); ++i ) {
				size += (( i > 0 ) ? joiner.size() : 0 ) + slices[i].size() ;
			}
			std::string result( size, 0 ) ;
			size = 0 ;
			for( std::size_t i = 0; i < slices.size(); ++i ) {
				std::copy( slices[i].begin(), slices[i].end(), result.begin() + size ) ;
				size += slices[i].size() ;
				std::copy( joiner.begin(), joiner.end(), result.begin() + size ) ;
				size += joiner.size() ;
			}
			return result ;
		}
		
		std::string operator+( slice const& left, slice const& right ) {
			std::string result( left.begin(), left.end() ) ;
			result += right ;
			return result ;
		}
		std::string operator+( slice const& left, std::string const& right ) {
			std::string result( left.begin(), left.end() ) ;
			result += right ;
			return result ;
		}
		std::string operator+( std::string const& left, slice const& right ) {
			std::string result( left ) ;
			result += right ;
			return result ;
		}
		std::string operator+( slice const& left, char const* right ) {
			std::string result( left ) ;
			result += std::string( right ) ;
			return result ;
		}
		std::string operator+( char const* left, slice const& right ) {
			std::string result( left ) ;
			result += std::string( right ) ;
			return result ;
		}
		
	}
}
