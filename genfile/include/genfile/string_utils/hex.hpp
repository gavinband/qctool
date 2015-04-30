
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_STRING_UTILS_HEX_HPP
#define GENFILE_STRING_UTILS_HEX_HPP

#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

namespace genfile {
	namespace string_utils {
		template< typename I, typename I2 >
		std::string to_hex( I i, I2 end_i ) {
			std::ostringstream o ;
			for( std::size_t count = 0; i < end_i; ++i, ++count ) {
				if( count % 4 == 0 )
					o << "|" ;
				o << std::hex << std::setw(2) << std::setfill('0') << static_cast<int> ( static_cast<unsigned char>( *i ) ) ;
			}
			return o.str() ;
		}

		std::string to_hex( std::string const& str ) ;
		std::string to_hex_char( std::string const& str ) ;
	}
}

#endif
