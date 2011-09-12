#include <iostream>
#include <string>
#include "../config.hpp"
#include "genfile/wildcard.hpp"
#include "statfile/statfile_utils.hpp"
#include "statfile/StatSink.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "statfile/BuiltInTypeStatSinkChain.hpp"

namespace statfile {
	std::auto_ptr< BuiltInTypeStatSinkChain > BuiltInTypeStatSinkChain::open( std::vector< genfile::wildcard::FilenameMatch > const& filenames ) {
		std::auto_ptr< BuiltInTypeStatSinkChain > chain( new BuiltInTypeStatSinkChain ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			chain->add_sink( BuiltInTypeStatSink::open( filenames[i].filename() )) ;
		}
		return chain ;
	}
	
	BuiltInTypeStatSinkChain::BuiltInTypeStatSinkChain(): m_current_sink(0) {}

	BuiltInTypeStatSinkChain::~BuiltInTypeStatSinkChain() {
		for( std::size_t i = 0; i < m_sinks.size(); ++i ) {
			delete m_sinks[i] ;
		}
	}

	void BuiltInTypeStatSinkChain::add_sink( std::auto_ptr< BuiltInTypeStatSink > sink ) {
		m_sinks.push_back( sink.release() ) ;
	}

	std::size_t BuiltInTypeStatSinkChain::index_of_current_sink() const {
		return m_current_sink ;
	}

	void BuiltInTypeStatSinkChain::add_column_impl( std::string const& name ) {
		base_t::add_column_impl( name ) ;
		for( std::size_t i = 0 ; i < m_sinks.size(); ++i ) {
			m_sinks[i]->add_column( name ) ;
		}
	}

	BuiltInTypeStatSinkChain::operator void*() const {
		if( m_current_sink < m_sinks.size() ) {
			return m_sinks[ m_current_sink ]->operator void*() ;
		}
		else {
			return 0 ;
		}
	}

	void BuiltInTypeStatSinkChain::write_value( int64_t const& value ) {
		if( *this ) {
			current_sink() << value ;
		}
	}

	void BuiltInTypeStatSinkChain::write_value( std::string const& value ) {
		if( *this ) {
			current_sink() << value ;
		}
	}

	void BuiltInTypeStatSinkChain::write_value( double const& value ) {
		if( *this ) {
			current_sink() << value ;
		}
	}
	
	void BuiltInTypeStatSinkChain::end_row() {
		if( *this ) {
			current_sink() << statfile::end_row() ;
		}
	}
	
	void BuiltInTypeStatSinkChain::move_to_next_sink() {
		assert( m_current_sink < m_sinks.size() ) ;
		++m_current_sink ;
	}

	std::size_t BuiltInTypeStatSinkChain::number_of_sinks() const { return m_sinks.size() ; }

	std::size_t BuiltInTypeStatSinkChain::get_current_sink() const {
		return m_current_sink ;
	}

	BuiltInTypeStatSink const& BuiltInTypeStatSinkChain::sink( std::size_t i ) const { return *m_sinks[ i ] ; }

	BuiltInTypeStatSink& BuiltInTypeStatSinkChain::current_sink() {
		assert( m_current_sink < m_sinks.size() ) ;
		return (*m_sinks[ m_current_sink ]) ;
	}
}

