#include <string>
#include <iostream>
#include "GToolException.hpp"

char GToolException::m_buffer[2000] ;

std::ostream& operator<<( std::ostream& out, GToolException const& exception ) {
    return out << GToolException::m_buffer ;
}


