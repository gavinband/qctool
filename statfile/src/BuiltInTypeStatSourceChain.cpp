
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include "statfile/BuiltInTypeStatSource.hpp"
#include "statfile/BuiltInTypeStatSourceChain.hpp"

namespace statfile {
	
	std::auto_ptr< BuiltInTypeStatSourceChain > BuiltInTypeStatSourceChain::open( std::vector< genfile::wildcard::FilenameMatch > const& filenames ) {
		std::auto_ptr< BuiltInTypeStatSourceChain > source( new BuiltInTypeStatSourceChain() ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			std::auto_ptr< BuiltInTypeStatSource > this_one = BuiltInTypeStatSource::open( filenames[i].filename() ) ;
			source->add_source( BuiltInTypeStatSource::open( filenames[i].filename() )) ;
		}
		return source ;
	}
	
	std::auto_ptr< statfile::BuiltInTypeStatSourceChain > BuiltInTypeStatSourceChain::open( std::vector< std::string > const& filenames ) {
		std::auto_ptr< BuiltInTypeStatSourceChain > source( new BuiltInTypeStatSourceChain() ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			source->add_source( BuiltInTypeStatSource::open( filenames[i] )) ;
		}
		return source ;
	}
	
	BuiltInTypeStatSourceChain::BuiltInTypeStatSourceChain(): m_current_source(0), m_moved_to_next_source_callback(0), m_number_of_rows(0) {}

	BuiltInTypeStatSourceChain::~BuiltInTypeStatSourceChain() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			delete m_sources[i] ;
		}
	}

	void BuiltInTypeStatSourceChain::add_source( std::auto_ptr< BuiltInTypeStatSource > source ) {
		if( m_sources.empty() ) {
			m_column_names = source->column_names() ;
			m_current_source = 0 ;
		}
		else if( source->column_names() != m_column_names ) {
			throw FileContainsRowsOfDifferentSizes() ;
		}
		m_sources.push_back( source.release() ) ;

		m_number_of_rows += m_sources.back()->number_of_rows() ;
	}

	void BuiltInTypeStatSourceChain::reset_to_start() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			m_sources[i]->reset_to_start() ;
		}
		m_current_source = 0 ;
		BuiltInTypeStatSource::reset_to_start() ;
	}

	std::size_t BuiltInTypeStatSourceChain::number_of_columns() const { return m_column_names.size() ; }
	std::vector< std::string > BuiltInTypeStatSourceChain::column_names() const { return m_column_names ; }
	std::string const& BuiltInTypeStatSourceChain::column_name( std::size_t i ) const { assert( i < m_column_names.size()) ; return m_column_names[i] ; }
	std::size_t BuiltInTypeStatSourceChain::number_of_rows() const {
		return m_number_of_rows ;
	}

	std::size_t BuiltInTypeStatSourceChain::number_of_sources() const { return m_sources.size() ; }

	unsigned int BuiltInTypeStatSourceChain::number_of_rows_in_source( std::size_t source_index ) const {
		assert( source_index < m_sources.size() ) ;
		return m_sources[ source_index ]->number_of_rows() ;
	}

	std::vector< unsigned int > BuiltInTypeStatSourceChain::get_source_row_counts() const {
		std::vector< unsigned int > result( number_of_sources() ) ;
		for( std::size_t i = 0; i < number_of_sources(); ++i ) {
			result[i] = number_of_rows_in_source( i ) ;
		}
		return result ;
	}

	void BuiltInTypeStatSourceChain::read_value( int32_t& value ) {
		move_to_next_nonempty_source_if_necessary() ;
		current_source() >> value ;
	}

	void BuiltInTypeStatSourceChain::read_value( uint32_t& value ) {
		move_to_next_nonempty_source_if_necessary() ;
		current_source() >> value ;
	}

	void BuiltInTypeStatSourceChain::read_value( std::string& value ) {
		move_to_next_nonempty_source_if_necessary() ;
		current_source() >> value ;
	}

	void BuiltInTypeStatSourceChain::read_value( double& value ) {
		move_to_next_nonempty_source_if_necessary() ;
		current_source() >> value ;
	}
	
	void BuiltInTypeStatSourceChain::ignore_value() {
		move_to_next_nonempty_source_if_necessary() ;
		current_source() >> IgnoreSome() ;
	}

	void BuiltInTypeStatSourceChain::ignore_all() {
		move_to_next_nonempty_source_if_necessary() ;
		current_source() >> IgnoreAll() ;
	}

	void BuiltInTypeStatSourceChain::end_row() {
		current_source() >> statfile::end_row() ;
	}

	void BuiltInTypeStatSourceChain::set_moved_to_next_source_callback( moved_to_next_source_callback_t callback ) { m_moved_to_next_source_callback = callback ; }

	BuiltInTypeStatSource& BuiltInTypeStatSourceChain::current_source() {
		assert( m_current_source < m_sources.size() ) ;
		return (*m_sources[ m_current_source ]) ;
	}

	void BuiltInTypeStatSourceChain::move_to_next_source() {
		++m_current_source ;
		if( m_moved_to_next_source_callback ) {
			m_moved_to_next_source_callback( m_current_source ) ;
		}
	}

	void BuiltInTypeStatSourceChain::move_to_next_nonempty_source_if_necessary() {
		while(( m_current_source < m_sources.size())
			&& (m_sources[ m_current_source ]->number_of_rows_read() >= m_sources[m_current_source]->number_of_rows())) {
			move_to_next_source() ;
		}
	}

}
