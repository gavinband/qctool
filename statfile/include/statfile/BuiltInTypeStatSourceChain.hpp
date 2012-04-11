
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef BUILTINTYPESTATSOURCECHAIN_HPP
#define BUILTINTYPESTATSOURCECHAIN_HPP

#include <iostream>
#include <string>
#include "../../config.hpp"
#if HAVE_BOOST_FUNCTION
#include <boost/function.hpp>
#endif
#include "genfile/wildcard.hpp"
#include "statfile/statfile_utils.hpp"
#include "statfile/StatSource.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"

namespace statfile {
	struct FileContainsRowsOfDifferentSizes: public FileStructureInvalidError { char const* what() const throw() { return "FileContainsRowsOfDifferentSizes" ; } } ;
	
	// class BuiltInTypeStatSourceChain represnets a BuiltInTypeStatSource
	// which gets it data sequentially from a collection of other StatSources
	class BuiltInTypeStatSourceChain: public ColumnNamingStatSource< BuiltInTypeStatSource >
	{
	public:
		typedef std::auto_ptr< statfile::BuiltInTypeStatSourceChain > UniquePtr ;
		static UniquePtr open( std::vector< genfile::wildcard::FilenameMatch > const& matches ) ;
		static UniquePtr open( std::vector< std::string > const& filenames ) ;
		
	public:
		
		BuiltInTypeStatSourceChain() ;
		~BuiltInTypeStatSourceChain() ;

		void add_source( std::auto_ptr< BuiltInTypeStatSource > source ) ;
		void reset_to_start() ;
		std::size_t number_of_columns() const ;
		std::vector< std::string > column_names() const ;
		std::string const& column_name( std::size_t i ) ;
		std::size_t number_of_rows() const ;
		std::size_t number_of_sources() const ;
		unsigned int number_of_rows_in_source( std::size_t source_index ) const ;
		std::vector< unsigned int > get_source_row_counts() const ;
		std::string get_descriptive_text() const { return "" ; }
		
	protected:
		void read_value( int32_t& value ) ;
		void read_value( uint32_t& value ) ;
		void read_value( std::string& value ) ;
		void read_value( double& value ) ;
		void ignore_value() ;
		void ignore_all() ;
		void end_row() ;

	protected:

	#if HAVE_BOOST_FUNCTION
		typedef boost::function< void( int index ) > moved_to_next_source_callback_t ;
	#else
		typedef void( *moved_to_next_source_callback_t )( std::size_t ) ;
	#endif

	public:
		void set_moved_to_next_source_callback( moved_to_next_source_callback_t callback ) ;

	private:

		BuiltInTypeStatSource& current_source() ;
		void move_to_next_source() ;
		void move_to_next_nonempty_source_if_necessary() ;

		std::vector< BuiltInTypeStatSource* > m_sources ;
		std::size_t m_current_source ;
		std::vector< std::string > m_column_names ;
		moved_to_next_source_callback_t m_moved_to_next_source_callback ;
		std::size_t m_number_of_rows ;
	} ;

}

#endif