
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_STRING_UTILS_SLICE_HPP
#define GENFILE_STRING_UTILS_SLICE_HPP

#include <vector>
#include <string>
#include <cassert>

namespace genfile {
	namespace string_utils {
		// class slice provides a (perhaps empty) view of a (possibly) larger string.
		// Assumption: the string object is managed elsewhere and outlives the view.
		// Invariant: any constructed slice is a valid view into a possibly
		// larger string.
		struct slice
		{
			// A view of the whole string.
			slice( std::string const& ) ;
			// A view of the whole string.
			slice( std::string const&, std::size_t start, std::size_t end ) ;
			slice( slice const&, std::size_t start, std::size_t end ) ;
			slice( slice const& other ) ;
			// slice& operator=( slice const& other ) ;

			char const& operator[]( std::size_t pos ) const { return (*m_string)[ m_start + pos ] ; }

			operator std::string() const { return m_string->substr( m_start, m_end - m_start ) ; }

			std::size_t size() const { return m_end - m_start ; }
			bool empty() const { return m_end == m_start ; }
			std::size_t get_start() const { return m_start ; }
			std::size_t get_end() const { return m_end ; }
			
			std::size_t find( char c, std::size_t pos = 0 ) const ;
			std::size_t find_first_of( std::string const& chars, std::size_t pos = 0 ) const ;

			std::vector< slice > split( std::string const& split_chars ) const ;
			void split( std::string const& split_chars, std::vector< slice >* result ) const ;
			
			bool operator==( std::string const& other ) const ;
			bool operator!=( std::string const& other ) const ;
			bool operator==( slice const& other ) const ;
		private:
			std::string const* m_string ;
			std::size_t m_start, m_end ;

		private:
			std::size_t find_first_of( char* membership_array, std::size_t pos = 0 ) const ;
		} ;
		
	}
}

#endif
