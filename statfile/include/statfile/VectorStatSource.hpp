
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef VECTOR_STAT_SOURCE_HPP
#define VECTOR_STAT_SOURCE_HPP

#include <vector>
#include "statfile/StatSource.hpp"

namespace statfile {
	template<
		typename T
	>
	struct VectorStatSource: public StatSource< std::size_t, T >
	{
		VectorStatSource( std::vector< T > const& data )
			: m_data( data )
		{
			build_column_names() ;
		}
	
		public:
			// Return the number of rows herein represented
			std::size_t number_of_rows() const { return m_data.size() ; }
			// Return the number of columns herein represented
			std::size_t number_of_columns() const { return m_column_names.size() ; }
			// Return the names of the columns
			std::vector< std::string > column_names() const {
				return m_column_names ;
			}

		void read_value( std::size_t& value ) {
			value = this->number_of_rows_read() ;
		}

		void read_value( T& value ) {
			assert( *this ) ;
			value = m_data[ this->number_of_rows_read() ] ;
		}

		protected:

			// In addition to the implementations of read( T& ) which are needed
			// for each type T != empty in the template list, derived classes must supply
			// the following two functions.
			void ignore_value() {} ;
			void ignore_all() {} ;
	
		private:
		
			void build_column_names() {
				m_column_names.push_back( "index" ) ;
				m_column_names.push_back( "value0" ) ;
			}
		
			std::vector< T > const m_data ;
			std::vector< std::string > m_column_names ;
	} ;
}

#endif
