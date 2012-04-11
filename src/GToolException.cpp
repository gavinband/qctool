
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <iostream>
#include "GToolException.hpp"

char GToolException::m_buffer[2000] ;

std::ostream& operator<<( std::ostream& out, GToolException const& exception ) {
    return out << GToolException::m_buffer ;
}


