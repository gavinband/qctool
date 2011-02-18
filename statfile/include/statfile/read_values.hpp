#ifndef STATFILE_READ_VALUES_HPP
#define STATFILE_READ_VALUES_HPP

#include <map>
#include <vector>
#include <string>
#include <boost/tuple/tuple.hpp>
#include "statfile/StatSource.hpp"

namespace statfile {
	namespace impl {
		std::map< std::size_t, std::size_t > get_indices_of_columns(
 			std::vector< std::string > const& source_column_names,
			std::string const& column_names
		) ;
		
		// To handle setting a value in a tuple, whose length we don't know, we make two functions.
		// The first takes a pointer to a character array of length 1 when the element is too large and 0 otherwise (and asserts).
		// The second takes a pointer to a character array of length 1 when the element is not too large.
		// Because char(*)[0] is not a valid type, only one of these can be instantiated.
		// This technique is an example of "SFINAE".
		template< int i, typename StatSource, typename ValueTuple >
		void read_value_into_tuple_elt( StatSource& source, ValueTuple const& values, char(*)[ i >= boost::tuples::length< ValueTuple >::value ] = 0 ) {
			assert(0) ; // not enough elements in tuple, we should never get here.
		}

		template< int i, typename StatSource, typename ValueTuple >
		void read_value_into_tuple_elt( StatSource& source, ValueTuple const& values, char(*)[ i < boost::tuples::length< ValueTuple >::value ] = 0 ) {
			source >> boost::tuples::get<i>( values ) ;
		}
	}
	
	template< typename StatSource, typename ValueTuple >
	bool read_values( StatSource& source, std::string const& column_names, ValueTuple const& values )
	// read_values is a convenience function which allows you to get data out of some columns
	// of a stat source using only the column names.  Column order issues are automatically
	// dealt with.
	{
		std::size_t const N = boost::tuples::length< ValueTuple >::value ;

		typedef std::map< std::size_t, std::size_t > IndexMap ;
		IndexMap column_indices ;
		try {
			 column_indices = impl::get_indices_of_columns( source.column_names(), column_names ) ;
		}
		catch( genfile::KeyNotFoundError const& e ) {
			throw genfile::KeyNotFoundError( "columns of " + source.get_source_spec(), e.key() ) ;
		}
		assert( column_indices.size() == N ) ;

		for( IndexMap::const_iterator i = column_indices.begin(); i != column_indices.end(); ++i ) {
			// get< j > needs a constant expression, so convert it to one here.
			source >> ignore( i->first - source.current_column() ) ;
			switch( i->second ) {
				case 0: impl::read_value_into_tuple_elt<0>( source, values ) ; break ;
				case 1: impl::read_value_into_tuple_elt<1>( source, values ) ; break ;
				case 2: impl::read_value_into_tuple_elt<2>( source, values ) ; break ;
				case 3: impl::read_value_into_tuple_elt<3>( source, values ) ; break ;
				case 4: impl::read_value_into_tuple_elt<4>( source, values ) ; break ;
				case 5: impl::read_value_into_tuple_elt<5>( source, values ) ; break ;
				case 6: impl::read_value_into_tuple_elt<6>( source, values ) ; break ;
				case 7: impl::read_value_into_tuple_elt<7>( source, values ) ; break ;
				case 8: impl::read_value_into_tuple_elt<8>( source, values ) ; break ;
				case 9: impl::read_value_into_tuple_elt<9>( source, values ) ; break ;
				default: assert(0) ; break ;
			}
			if( !source ) {
				return false ;
			}
		}
		return source ;
	}
}

#endif

