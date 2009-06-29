#ifndef __GTOOL_ObjectSOURCE_HPP__
#define __GTOOL_ObjectSOURCE_HPP__

#include <cassert>

// Base class for classes from which Objects can be read.
template< typename Object >
struct ObjectSource
{
public:
	virtual ~ObjectSource() {};
	
	virtual void read( Object& ) = 0 ;
	virtual bool check_if_empty() = 0 ;
	virtual operator bool() = 0 ;
} ;

template< typename Object >
ObjectSource< Object >& operator>>( ObjectSource< Object >& source, Object& object ) {
	assert( !source.check_if_empty() ) ;
	source.read( object ) ;
	return source ;
}

template< typename Object >
struct NullObjectSource: public ObjectSource< Object >
{
	void read( Object& ) { assert( 0 ) ;}
	bool check_if_empty() { return true ;}
	operator bool() { return false ;}
} ;


#endif
