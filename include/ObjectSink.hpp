#ifndef __GTOOL_ObjectSink_HPP__
#define __GTOOL_ObjectSink_HPP__

#include <cassert>

// Base class for classes to which GenRows can be written.
template< typename Object >
struct ObjectSink
{
public:
	
	virtual ~ObjectSink() {};
	
	virtual void write( Object const& ) = 0;
	virtual bool check_if_full() = 0;
} ;

template< typename Object >
ObjectSink<Object> & operator<<( ObjectSink< Object >& sink, Object const& object ) {
	assert( !sink.check_if_full()) ;
	sink.write( object ) ;
	return sink ;
}

template< typename Object >
struct NullObjectSink: public ObjectSink< Object >
{
	void write( Object const& ) {}
	bool check_if_full() { return false ;}
} ;

#endif