#ifndef VECTOR_OF_VECTORS_STAT_SOURCE_HPP
#define VECTOR_OF_VECTORS_STAT_SOURCE_HPP

#include <vector>
#include "StatSource.hpp"

template<
	typename T
>
struct VectorOfVectorsStatSource: public StatSource< std::size_t, T >
{
	VectorStatSource( std::vector< std::vector< T > > const& data )
		: m_data( data )
	{
		build_column_names() ;
		assert( verify_row_sizes() ) ;
	}
	
	public:
		// Return the number of rows herein represented
		std::size_t number_of_rows() const { return m_data.size() ; }
		// Return the number of columns herein represented
		std::size_t number_of_columns() const {
			return m_column_names.size() ;
		}
		// Return the names of the columns
		std::vector< std::string > column_names() const {
			return m_column_names ;
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
		
		bool verify_row_sizes() const {
			bool result = true ;

			if( !m_data.empty() ) {
				std::size_t row_size = m_data[0].size() ;
				for( std::size_t i = 1; i < m_data.size(); ++i ) {
					if( m_data[i].size() != row_size ) {
						result = false ;
					}
				}
			}
			return result ;
		}
		
		void build_column_names() {
			m_column_names.push_back( "index" ) ;
			if( m_data.size() > 0 ) {
				for( std::size_t i = 0; i < m_data[0].size(); ++i ) {
					m_column_names.push_back( "value" + to_string( i )) ;
				}
			}
		}

		std::string to_string( std::size_t i ) {
			std::ostringstream ostr ;
			ostr << i ;
			return ostr.str() ;
		}
		
		std::vector< std::vector< T > > const& m_data ;
		std::vector< std::string > m_column_names ;
} ;

#endif
