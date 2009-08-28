#ifndef PROGRAM_FLOW_HPP
#define PROGRAM_FLOW_HPP

#include <exception>

struct HaltProgramWithReturnCode: public std::exception {
	HaltProgramWithReturnCode( int return_code = -1 ) : m_return_code( return_code ) {}
	~HaltProgramWithReturnCode() throw() {}
	char const* what() const throw() { return "HaltProgramWithReturnCode" ; }
	int return_code() const { return m_return_code ; }
private:
	int m_return_code ;
} ;

#endif
