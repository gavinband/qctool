#include <string>
#include <ctime>
#include "appcontext/get_current_time_as_string.hpp"

namespace appcontext {
	std::string get_current_time_as_string() {
		time_t rawtime ;
		struct tm * timeinfo ;
		char buffer[30] ;
		std::time( &rawtime ) ;
		timeinfo = std::localtime( &rawtime ) ;
		std::strftime( buffer, 80, "%Y-%m-%d %H:%M:%S", timeinfo ) ;
		return std::string( buffer ) ;
	}
}
