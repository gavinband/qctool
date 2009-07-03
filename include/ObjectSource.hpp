#ifndef __GTOOL_ObjectSOURCE_HPP__
#define __GTOOL_ObjectSOURCE_HPP__

#include <cassert>

// Base class for classes from which Objects can be read.
template< typename Object >
struct ObjectSource
{
public:
	virtual ~ObjectSource() {};
	
	virtual ObjectSource& read( Object& ) = 0 ;
	virtual operator bool() = 0 ;
	virtual bool fail() const = 0 ;
} ;

template< typename Object >
ObjectSource< Object >& operator>>( ObjectSource< Object >& source, Object& object ) {
	source.read( object ) ;
	return source ;
}

template< typename Object >
struct NullObjectSource: public ObjectSource< Object >
{
	NullObjectSource& read( Object& ) { return *this ; }
	operator bool() { return false ;}
	bool fail() const { return false ; }
} ;


#endif
