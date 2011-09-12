#ifndef BUILTINTYPESTATSINKCHAIN_HPP
#define BUILTINTYPESTATSINKCHAIN_HPP

#include <iostream>
#include <string>
#include "../config.hpp"
#include "genfile/wildcard.hpp"
#include "statfile/statfile_utils.hpp"
#include "statfile/StatSink.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"

namespace statfile {
	// class BuiltInTypeStatSinkChain represnets a BuiltInTypeStatSink
	// which gets it data sequentially from a collection of other StatSources
	class BuiltInTypeStatSinkChain: public ColumnNamingStatSink< BuiltInTypeStatSink >
	{
		typedef ColumnNamingStatSink< BuiltInTypeStatSink > base_t ;
	public:
		static std::auto_ptr< BuiltInTypeStatSinkChain > open( std::vector< genfile::wildcard::FilenameMatch > const& ) ;
	public:
		
		BuiltInTypeStatSinkChain() ;
		~BuiltInTypeStatSinkChain() ;
		void add_sink( std::auto_ptr< BuiltInTypeStatSink > sink ) ;
		std::size_t index_of_current_sink() const ;
		void add_column_impl( std::string const& name ) ;
		operator void*() const ;
		void write_value( int32_t const& value ) ;
		void write_value( int64_t const& value ) ;
		void write_value( uint32_t const& value ) ;
		void write_value( std::string const& value ) ;
		void write_value( double const& value ) ;
		void end_row() ;
		void move_to_next_sink() ;
		std::size_t number_of_sinks() const ;
		BuiltInTypeStatSink const& sink( std::size_t i ) const ;
		std::size_t get_current_sink() const ;
	private:

		BuiltInTypeStatSink& current_sink() ;
		std::vector< BuiltInTypeStatSink* > m_sinks ;
		std::size_t m_current_sink ;
		std::vector< std::string > m_column_names ;
	} ;
}

#endif