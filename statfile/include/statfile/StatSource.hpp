#ifndef PERSNPDATASOURCE_HPP
#define PERSNPDATASOURCE_HPP

#include <iostream>
#include <string>
#include <vector>
#include <boost/function.hpp>
#include "statfile/statfile_utils.hpp"
#include "statfile/types.hpp"
#include "statfile/TypeReaderBase.hpp"

namespace statfile {
	struct StatSourceError: public StatError { char const* what() const throw() { return "StatSourceError" ; } } ;
	struct FileIsInvalidError: public StatSourceError { char const* what() const throw() { return "FileIsInvalidError" ; } } ;
	struct FileHasTrailingDataAfterLastSNP: public FileIsInvalidError { char const* what() const throw() { return "FileHasTrailingDataAfterLastSNP" ; } } ;
	struct FileContainsSNPsOfDifferentSizes: public FileIsInvalidError { char const* what() const throw() { return "FileContainsSNPsOfDifferentSizes" ; } } ;


	// Base class for classes which provide SNP assay data, one snp at a time.
	// After the class is constructed, the intention is that
	// 1. the number_of_numerical_entries_per_snp() and total_number_of_snps() functions return information
	// reflecting the data in the file or source, and
	// 2. The data can be read using read() (or get_identifying_data() and read_statfile())
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
		StatSource(): m_current_column(0u), m_number_read(0u) {} ;
		virtual ~StatSource() {} ;

	private:
		// We forbid copying & assignment.
		StatSource( StatSource const& other ) ;
		StatSource& operator=( StatSource const& other ) ;

	public:

		template< typename T >
		StatSource& operator>>( T& value ) {
			assert( m_current_column < number_of_columns() ) ;
			this->read_value( value ) ;
			if( *this ) {
				move_to_next_column() ;
			}
			return *this ;
		}

		StatSource& operator>>( IgnoreOne const& ) {
			assert( m_current_column < number_of_columns() ) ;
			this->ignore_value() ;
			if( *this ) {
				move_to_next_column() ;
			}
			return *this ;
		}

		StatSource& operator>>( IgnoreAll const& ) {
			assert( m_current_column < number_of_columns() ) ;
			this->ignore_all() ;
			if( *this ) {
				move_to_next_row() ;
			}
			return *this ;
		}


	public:
		// Return nonzero pointer iff there have been no errors up to now.
		virtual operator void*() const = 0 ;
		// Return the number of rows herein represented
		virtual std::size_t total_number_of_rows() const = 0 ;
		// Return the number of rows read so far
		std::size_t number_of_rows_read() const { return m_number_read ; }
		// Return the number of columns herein represented
		virtual std::size_t number_of_columns() const = 0;
		// Return the names of the columns
		virtual std::vector< std::string > const& column_names() const = 0 ;

	protected:

		// In addition to the implementations of read( T& ) which are needed
		// for each type T != empty in the template list, derived classes must supply
		// the following two functions.
		virtual void ignore_value() = 0 ;
		virtual void ignore_all() = 0 ;

		// Return the current column.
		std::size_t current_column() { return m_current_column ; }

	private:
		void move_to_next_column() {
			if( ++m_current_column == number_of_columns() ) {
				move_to_next_row() ;
			}
		}
		
		void move_to_next_row() {
			move_to_next_row_impl() ;
			m_current_column = 0 ;
			++m_number_read ;
		}
		
		virtual void move_to_next_row_impl() {} ;
		
		std::size_t m_current_column ;
		std::size_t m_number_read ; // keep track of #of snps read from the source
		
	} ;

	template<
		typename StatSourceT
	>
	struct ColumnNamingStatSource: public virtual StatSourceT
	{
		std::size_t number_of_columns() const { return m_column_names.size() ; }
		std::vector< std::string > const& column_names() const { return m_column_names ; }
	protected:
		void add_column( std::string const& name ) {
			m_column_names.push_back( name ) ;
		}
	private:		
		std::vector< std::string > m_column_names ;
	} ;
	
	typedef StatSource< int32_t, uint32_t, std::string, double > BuiltInTypeStatSource ;
}



#endif