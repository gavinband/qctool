#include <exception>
#define GENFILE_DECLARE_EXCEPTION( name ) 							\
struct name: public std::exception 									\
{																	\
	char const* what() const throw() { return "name" ; }			\
} ;																	\
// end of GENFILE_DECLARE_EXCEPTION

#define GENFILE_DECLARE_SINGLE_ARG_EXCEPTION( name, type, arg ) 	\
struct name: public std::exception 									\
{																	\
	name( type const& arg ): m_arg( arg ) {}						\
	~name() throw() {}												\
	char const* what() const throw() { return "name" ; }			\
	type const& arg() const { return m_arg ; }						\
private																\
	type const m_arg ;												\
} ;																	\
// end of GENFILE_DECLARE_EXCEPTION_WITH_ARG
