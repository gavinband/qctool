
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef STATFILE_FROM_STRING_HPP
#define STATFILE_FROM_STRING_HPP

#include <string>

namespace statfile {
	template< typename T >
	class FromString
	{
	public:
		FromString( T& t ):
			m_t( t )
		{}
		
		~FromString() {
			m_t = T( m_string ) ;
		}
		
		operator std::string& () { return m_string ;}
		
	private:
		T& m_t ;
		std::string m_string ;
	}  ;
	
	template< typename T >
	FromString< T > from_string( T& t ) {
		return FromString< T >( t ) ;
	}
}

#endif
