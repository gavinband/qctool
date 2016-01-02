
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_STRING_UTILS_SLICE_HPP
#define GENFILE_STRING_UTILS_SLICE_HPP

#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <boost/function.hpp>

namespace genfile {
	namespace string_utils {
		// class slice provides a (perhaps empty) view of a (possibly) larger string.
		// Assumption: the string object is managed elsewhere and outlives the view.
		// Invariant: any constructed slice is a valid view into a possibly
		// larger string.
		struct slice
		{
		public:
			typedef char const* const_iterator ;

		public:
			// A view of a C-style string; view does not include terminating null.
			slice( char const* c_string ) ;
			// A view of a subset of a c-style string.
			slice( char const* c_string, std::size_t start, std::size_t end ) ;
			// A view of a character buffer
			slice( char const* data, char const* end_of_data ) ;
			// A view of a subset of a character buffer
			slice( char const* data, char const* end_of_data, std::size_t start, std::size_t end ) ;
			// A view of a string.
			slice( std::string const& ) ;
			// A view of a subset of a string.
			slice( std::string const&, std::size_t start, std::size_t end ) ;
			slice( slice const&, std::size_t start, std::size_t end ) ;
			slice( slice const& other ) ;
			// slice& operator=( slice const& other ) ;

			char const& operator[]( std::size_t pos ) const { return *(m_data + m_start + pos) ; }
			operator std::string() const { return std::string( m_data + m_start, m_data + m_end ) ; }

			std::size_t size() const { return m_end - m_start ; }
			bool empty() const { return m_end == m_start ; }
			std::size_t get_start() const { return m_start ; }
			std::size_t get_end() const { return m_end ; }
			
			std::size_t find( char c, std::size_t pos = 0 ) const ;
			std::size_t find_first_of( std::string const& chars, std::size_t pos = 0 ) const ;
			std::size_t find_first_not_of( std::string const& chars, std::size_t pos = 0 ) const ;
			std::size_t find_last_of( std::string const& chars, std::size_t pos = std::string::npos ) const ;
			std::size_t find_last_not_of( std::string const& chars, std::size_t pos = std::string::npos ) const ;

			slice strip( std::string const& strip_chars ) const ;

			std::vector< slice > split( std::string const& split_chars ) const ;
			void split( std::string const& split_chars, std::vector< slice >* result ) const ;
			void split( std::string const& split_chars, boost::function< void( slice ) > ) const ;
			
			slice substr( std::size_t start, std::size_t end ) const ;
			
			bool operator==( std::string const& other ) const ;
			bool operator!=( std::string const& other ) const ;
			bool operator==( char const* other ) const ;
			bool operator!=( char const* other ) const ;
			// bool operator==( slice const& other ) const ;
			
			const_iterator begin() const ;
			const_iterator end() const ;

			friend bool operator<( slice const& left, slice const& right ) ;
			friend bool operator>( slice const& left, slice const& right ) ;
			friend bool operator<=( slice const& left, slice const& right ) ;
			friend bool operator>=( slice const& left, slice const& right ) ;
			friend bool operator==( slice const& left, slice const& right ) ;
			friend bool operator!=( slice const& left, slice const& right ) ;
			friend std::ostream& operator<<( std::ostream& o, slice const& s ) ;
			
		private:
			char const* const m_data ;
			char const* const m_end_of_data ;
			std::size_t m_start, m_end ;

		private:
			std::size_t find_first_of( char* membership_array, std::size_t pos = 0 ) const ;
			std::size_t find_last_of( char* membership_array, std::size_t pos = std::string::npos ) const ;
		} ;
	
		std::string join( std::vector< slice > const& slices, std::string const& joiner ) ;

		std::string operator+( slice const& left, slice const& right ) ;
		std::string operator+( slice const& left, std::string const& right ) ;
		std::string operator+( std::string const& left, slice const& right ) ;
		std::string operator+( slice const& left, char const* right ) ;
		std::string operator+( char const* left, slice const& right ) ;
	}
}

#endif
