#ifndef __GTOOL_CHAININGFILEOBJECTSOURCE_HPP
#define __GTOOL_CHAININGFILEOBJECTSOURCE_HPP


#include <vector>
#include "ObjectSource.hpp"

template< typename Object >
struct ChainingFileObjectSource: public ObjectSource< Object >
{
	ChainingFileObjectSource()
		: m_current_source(0)
	{}

	~ChainingFileObjectSource() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			delete m_sources[i] ;
		}
	}

	void add_source( std::auto_ptr< ObjectSource< Object > >& source ) {
		m_sources.push_back( source.release() ) ;
	}

	ChainingFileObjectSource& read( Object& row ) {
		while(( m_current_source < m_sources.size()) && (! (m_sources[ m_current_source ]->read( row )))) {
			++m_current_source ;
		}
		return *this ;
	}
	
	operator bool() {
		for( ; m_current_source < m_sources.size(); ++m_current_source ) {
			if( *(m_sources[ m_current_source ])) {
				return true ;
			}
		}
		
		return false ;
	}

	bool fail() const {
		if( m_current_source < m_sources.size() ) {
			return m_sources[ m_current_source ]->fail() ;
		}
		else {
			return false ;
		}
	}

	std::size_t current_source() const { return m_current_source ; }
	std::size_t number_of_sources() const { return m_sources.size() ; }

private:
	
	void move_to_next_nonempty_source() {
		while( m_current_source < m_sources.size() && m_sources[m_current_source]->check_if_empty() ) {
			++m_current_source ;
		}
	}

	std::size_t m_current_source ;
	std::vector< ObjectSource< Object >* > m_sources ;
} ;

#endif
