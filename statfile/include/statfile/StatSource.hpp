
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef STATFILE_STATSOURCE_HPP
#define STATFILE_STATSOURCE_HPP

#include <iostream>
#include <string>
#include <vector>
#include <boost/function.hpp>
#include "statfile/statfile_utils.hpp"
#include "statfile/types.hpp"
#include "statfile/TypeReaderBase.hpp"

namespace statfile {
	struct StatSourceError: public StatError { char const* what() const throw() { return "StatSourceError" ; } } ;
	struct FileHasTrailingDataAfterLastSNP: public StatSourceError { char const* what() const throw() { return "FileHasTrailingDataAfterLastSNP" ; } } ;
	struct FileContainsSNPsOfDifferentSizes: public StatSourceError { char const* what() const throw() { return "FileContainsSNPsOfDifferentSizes" ; } } ;
	struct ReadPastEndError: public StatSourceError { char const* what() const throw() { return "ReadPastEndError" ; } } ;


	// Base class for classes which provide data, nominally organised into rows and columns.
	// Data can be read using operator>>, which reads the next entry from the current row.
	// You must call source >> end_row() at the end of each row -- this helps make sure the user
	// is in sync with the source.
	//
	// Methods are also provided to return the current number of rows and columns, and the number
	// of rows read from the source.
	//
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
	class StatSource: public TypeReaderBase< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 >
	{
	public:
		StatSource(): m_current_column(0u), m_number_read(0u), m_source_exhausted( false ) {} ;
		virtual ~StatSource() {} ;

		virtual void reset_to_start() {
			m_current_column = 0u ;
			m_number_read = 0u ;
			m_source_exhausted = false ;
		}

	private:
		// We forbid copying & assignment.
		StatSource( StatSource const& other ) ;
		StatSource& operator=( StatSource const& other ) ;

	public:

		template< typename T >
		StatSource& operator>>( T& value ) {
			if( m_number_read == number_of_rows() ) {
				m_source_exhausted = true ;
			}
			else {
				assert( m_current_column < number_of_columns()) ;
				this->read_value( value ) ;
				move_to_next_column() ;
			}
			return *this ;
		}

		StatSource& operator>>( IgnoreSome const& ignore_some ) {
			if( m_number_read == number_of_rows() ) {
				m_source_exhausted = true ;
			}
			else {
				std::size_t target_column = m_current_column + ignore_some.number_to_ignore() ;
				assert( target_column < number_of_columns() ) ;
				for( std::size_t i = 0; i < ignore_some.number_to_ignore() ; ++i ) {
					this->ignore_value() ;
					move_to_next_column() ;
				}
			}
			return *this ;
		}

		StatSource& operator>>( IgnoreAll const& ) {
			if( m_number_read == number_of_rows() ) {
				m_source_exhausted = true ;
			}
			else {
				assert( m_current_column <= number_of_columns() ) ;
				this->ignore_all() ;
				move_to_next_row() ;
			}
			return *this ;
		}

		// TODO: make use of EndRow mandatory.
		StatSource& operator>>( EndRow const& ) {
			if( m_number_read == number_of_rows() ) {
				m_source_exhausted = true ;
			}
			else {
				assert( m_current_column == number_of_columns() ) ;
				end_row() ;
				move_to_next_row() ;
			}
			return *this ;
		}

	public:
		// Return nonzero pointer iff there have been no errors up to now.
		operator void const*() const { return m_source_exhausted ? 0 :  reinterpret_cast< void const* >( this ) ; }
		// Return the number of rows herein represented
		virtual std::size_t number_of_rows() const = 0 ;
		// Return the number of rows read so far
		std::size_t number_of_rows_read() const { return m_number_read ; }
		// Return the number of columns herein represented
		virtual std::size_t number_of_columns() const = 0;
		// Return the names of the columns
		virtual std::vector< std::string > column_names() const = 0 ;
		// Return the name of the ith column
		std::string name_of_column( std::size_t index ) const {
			std::vector< std::string > const& names = column_names() ;
			assert( index < names.size() ) ;
			return names[ index ] ;
		} ;
		// Return the index of the named column.  Index must lie in the range [ 0, number_of_columns() ).
		std::size_t index_of_column( std::string const& name ) const {
			std::vector< std::string > const& names = column_names() ;
			std::vector< std::string >::const_iterator where = std::find( names.begin(), names.end(), name ) ;
			if( where == names.end() ) {
				throw InvalidColumnNameError( name ) ;
			}
			return where - names.begin() ;
		}

		// Return true or false if the source has the given column.
		bool has_column( std::string const& name ) const {
			std::vector< std::string > const& names = column_names() ;
			return( std::find( names.begin(), names.end(), name ) != names.end() ) ;
		}

		// Return the current column.
		std::size_t current_column() { return m_current_column ; }
		virtual std::string get_source_spec() const { return "(unknown)" ; }

	protected:

		// In addition to the implementations of read( T& ) which are needed
		// for each type T != empty in the template list, derived classes must supply
		// the following two functions.
		virtual void ignore_value() = 0 ;
		virtual void ignore_all() = 0 ;
		virtual void end_row() {} ;

	private:
		void move_to_next_column() {
			assert( ++m_current_column <= number_of_columns() ) ;
		}

		void move_to_column() {
			assert( ++m_current_column <= number_of_columns() ) ;
		}
		
		void move_to_next_row() {
			move_to_next_row_impl() ;
			m_current_column = 0 ;
			++m_number_read ;
		}
		
		virtual void move_to_next_row_impl() {} ;
		
		std::size_t m_current_column ;
		std::size_t m_number_read ; // keep track of #of snps read from the source

		bool m_source_exhausted ;	
	} ;

	template< typename StatSourceT >
	struct ColumnNamingStatSourceBase: public virtual StatSourceT
	{
		virtual std::string const& column_name( std::size_t index ) const = 0 ;
	protected:
		virtual void add_column( std::string const& ) = 0 ;
	} ;

	template<
		typename StatSourceT
	>
	struct ColumnNamingStatSource: public ColumnNamingStatSourceBase< StatSourceT >
	{
		std::size_t number_of_columns() const { return m_column_names.size() ; }
		std::vector< std::string > column_names() const { return m_column_names ; }
		std::string const& column_name( std::size_t index ) const {
			assert( index < m_column_names.size() ) ;
			return m_column_names[ index ] ;
		}
	protected:
		void add_column( std::string const& name ) {
			m_column_names.push_back( name ) ;
		}
	private:		
		std::vector< std::string > m_column_names ;
	} ;
}



#endif

