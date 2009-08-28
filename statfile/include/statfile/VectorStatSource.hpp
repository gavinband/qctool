#ifndef VECTOR_STAT_SOURCE_HPP
#define VECTOR_STAT_SOURCE_HPP

#include <vector>
#include "StatSource.hpp"

template<
	typename T
>
struct VectorStatSource: public StatSource< std::size_t, T >
{
	VectorStatSource( std::vector< T > const& data )
		: m_data( data )
	{}
	
	public:
		// Return nonzero pointer iff there have been no errors up to now.
		operator void*() const { return (number_of_rows_read() < m_data.size()) ? this : 0 ; }
		// Return the number of rows herein represented
		std::size_t total_number_of_rows() const { return m_data.size() ; }
		// Return the number of columns herein represented
		std::size_t number_of_columns() const { return 2 ; }
		// Return the names of the columns
		std::vector< std::string > const& column_names() const {
			std::vector< std::string > result ;
			result.push_back( "index" ) ;
			result.push_back( "value" ) ;
			return result ;
		}

	void read_value( std::size_t& value ) {
		return number_of_rows_read() ;
	}

	void read_value( T& value ) {
		return m_data[ number_of_rows_read() ] ;
	}

	protected:

		// In addition to the implementations of read( T& ) which are needed
		// for each type T != empty in the template list, derived classes must supply
		// the following two functions.
		void ignore_value() {} ;
		void ignore_all() {} ;
	
	private:
		
		std::vector< T > const& m_data ;
} ;

#endif
