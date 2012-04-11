
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GEN_TOOLS_R_HPP
#define GEN_TOOLS_R_HPP

#include "../config.hpp"
#include <stdlib.h>

namespace R {
	int interpret_R_code( std::string const& R_code ) {
		int result = 0 ;
		if( system( NULL )) {
			std::string tmp_filename = tmpnam(0) ;
			std::ofstream R_command_file( tmp_filename.c_str() ) ;
			if( !R_command_file.is_open() ) {
				return -1 ;
			}
			R_command_file << R_code ;
			R_command_file.close() ;
			
			std::string command = "/usr/bin/R -f \"" + tmp_filename + "\" >> /dev/null" ;
			result = system( command.c_str() ) ;
		}
		
		return result ;
	}
}


#endif
