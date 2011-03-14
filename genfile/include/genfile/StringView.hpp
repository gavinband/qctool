#ifndef GENFILE_STRINGVIEW_HPP
#define GENFILE_STRINGVIEW_HPP

#include <vector>
#include <string>

namespace genfile {
	namespace stringview {
		// class StringView provides a (perhaps empty) view of a (possibly) larger string.
		// Assumption: the string object is managed elsewhere and outlives the view.
		// Invariant: any constructed StringView is a valid view into a possibly
		// larger string.
		struct StringView
		{
			// A view of the whole string.
			StringView( std::string const& ) ;
			// A view of the whole string.
			StringView( std::string const&, std::size_t start, std::size_t end ) ;
			StringView( StringView const&, std::size_t start, std::size_t end ) ;
			StringView( StringView const& other ) ;
			StringView& operator=( StringView const& other ) ;

			char operator[]( std::size_t pos ) const { return (*m_string)[ m_start + pos ] ; }
			char const& operator[]( std::size_t pos ) { return (*m_string)[ m_start + pos ] ; }

			operator std::string() const { return m_string->substr( m_start, m_end - m_start ) ; }

			std::size_t size() const { return m_end - m_start ; }
			std::size_t get_start() const { return m_start ; }
			std::size_t get_end() const { return m_end ; }
			
			std::size_t find( char c, std::size_t pos = 0 ) const ;
			std::size_t find_first_of( std::string const& chars, std::size_t pos = 0 ) const ;
			
			bool operator==( std::string const& other ) const ;
		private:
			std::string const* m_string ;
			std::size_t m_start, m_end ;
		} ;
		
		std::vector< StringView > split( StringView const& string_to_split, std::string const& split_chars ) ;
	}
}

#endif
