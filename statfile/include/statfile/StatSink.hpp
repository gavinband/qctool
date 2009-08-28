#ifndef STATSINK_HPP
#define STATSINK_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <stddef.h>
#include <boost/function.hpp>
#include "statfile/statfile_utils.hpp"
#include "statfile/TypeWriterBase.hpp"

namespace statfile {
	struct StatSinkError: public std::exception { char const* what() const throw() { return "StatSinkError" ; } } ;

	// Base class representing a sink of data which is stored in the following format.
	// Each (conceptual) row must have an index field (a 32-bit unsigned integer) followed by
	// a number of "identifying fields" (each a string), followed by the data itself.
	template<
		typename T1 = empty,
		typename T2 = empty,
		typename T3 = empty,
		typename T4 = empty,
		typename T5 = empty,
		typename T6 = empty,
		typename T7 = empty,
		typename T8 = empty,
		typename T9 = empty,
		typename T10 = empty
	>
	class StatSink: public TypeWriterBase< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 >
	{
	public:
		StatSink()
		 : m_state( e_HaveNotWrittenAnyData ),
		   m_current_column(0u),
		   m_number_of_rows_written(0u)
		{}
	private:
		// Forbid copying and assignment
		StatSink( StatSink const& other ) ;
		StatSink& operator=( StatSink const& other ) ;
	public:
		void add_columns( std::vector< std::string > const& names ) {
			for( std::size_t i = 0; i < names.size(); ++i ) {
				add_column( names[i] ) ;
			}
		} ;

		// Use this to add index and identifying members to each row.
		template< typename T >
		StatSink& operator<<( T const& value ) {
			assert( m_current_column < number_of_columns() ) ;
			write_value( value ) ;
			m_state = e_HaveWrittenSomeData ;
			if( *this ) {
				move_to_next_column() ;
			}
			return *this ;
		}

		StatSink& operator<<( EndRow const& ) {
			assert( m_current_column == number_of_columns()) ;
			move_to_next_row() ;
			return *this ;
		}

	public:
		// In addition to write( T const& ) which derived classes must supply for
		// each T != empty in the template parameter list, the following functions
		// must be supplied by each derived class.

		// Return a nonzero pointer iff there have been no errors so far.
		virtual operator void*() const = 0 ;
		// Column-related functions which must be supplied by derived classes.
		virtual void add_column( std::string const& name ) = 0 ;
		virtual std::size_t number_of_columns() const = 0 ;
		virtual std::vector< std::string > const& column_names() const = 0 ;

		// The next two functions return the tracked column and row numbers.
		std::size_t number_of_rows_written() const { return m_number_of_rows_written ; }
		std::size_t current_column() const { return m_current_column ; }

	protected:
		
		enum State {
			e_HaveNotWrittenAnyData = 0,
			e_HaveWrittenSomeData = 1
		} ;

		int state() const { return m_state ; }

	private:
		void move_to_next_column() {
			assert( m_current_column < number_of_columns() ) ;
			++m_current_column ;
		}

		void move_to_next_row() {
			move_to_next_row_impl() ;
			m_current_column = 0 ;
			++m_number_of_rows_written ;
		}
	
		virtual void move_to_next_row_impl() {} ;
	
		int m_state ;
		std::size_t m_current_column ;
		std::size_t m_number_of_rows_written ;
	} ;
	
	typedef StatSink< int, unsigned int, long int, long unsigned int, std::string, double > BuiltInTypeStatSink ;
	
	template<
		typename StatSinkT
	>
	struct ColumnNamingStatSink: public virtual StatSinkT
	{
		void add_column( std::string const& name ) {
			m_column_names.push_back( name ) ;
		} ;
		std::size_t number_of_columns() const { return m_column_names.size(); }
		std::vector< std::string > const& column_names() const { return m_column_names ; }
	private:
		std::vector< std::string > m_column_names ;
	} ;
	
	struct TrivialBuiltInTypeStatSink: public ColumnNamingStatSink< BuiltInTypeStatSink >
	{
		operator void*() const { return reinterpret_cast< void* > ( const_cast< TrivialBuiltInTypeStatSink* >( this )) ; }
		void write_value( int const& ) {}
		void write_value( unsigned int const& ) {}
		void write_value( long int const& ) {}
		void write_value( long unsigned int const& ) {}
		void write_value( double const& ) {}
		void write_value( std::string const& ) {}
	} ;
}

#endif