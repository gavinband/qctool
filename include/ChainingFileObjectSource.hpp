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

	bool check_if_empty() {
		find_next_nonempty_source() ;
		return m_current_source >= m_sources.size() ;
	}

	void read( Object& row ) {
		find_next_nonempty_source() ;
		assert( !(m_sources[ m_current_source ]->check_if_empty())) ;
		m_sources[ m_current_source ]->read( row ) ;
	}

	std::size_t current_source() const { return m_current_source ; }
	std::size_t number_of_sources() const { return m_sources.size() ; }

private:
	
	void find_next_nonempty_source() {
		while( m_current_source < m_sources.size() && m_sources[m_current_source]->check_if_empty() ) {
			++m_current_source ;
		}
	}

	std::size_t m_current_source ;
	std::vector< ObjectSource< Object >* > m_sources ;
} ;

#endif
