#ifndef GENFILE_UTILITY_HPP
#define GENFILE_UTILITY_HPP

#include <vector>
#include <algorithm>
#include <iterator>

namespace genfile {
	namespace utility {
		// function: select_entries
		// Return a vector having only the entries from the given vector whose index lies in the
		// second argument.
		template< typename V >
		std::vector< V > select_entries( std::vector< V > const& container, std::vector< std::size_t > const& indices ) {
			std::vector< V > result( indices.size() ) ;
			for( std::size_t i = 0; i < indices.size(); ++i ) {
				assert( indices[i] < container.size() ) ;
				result[i] = container[ indices[i] ] ;
			}
			return result ;
		}
		
		// function: select_entries
		// Return a vector having only those entries of the first argument for which the corresponding entry of the
		// second argument is true.
		template< typename V >
		std::vector< V > select_entries( std::vector< V > const& container, std::vector< bool > const& filter ) {
			assert( container.size() == filter.size() ) ;
			std::vector< V > result ;
			result.reserve( container.size() ) ;
			for( std::size_t i = 0; i < container.size(); ++i ) {
				if( filter[i] ) {
					result.push_back( container[i] ) ;
				}
			}
			return result ;
		}
		
		template< typename S >
		S intersect( S const& set1, S const& set2 ) {
			S result ;
			std::set_intersection(
				set1.begin(), set1.end(),
				set2.begin(), set2.end(),
				std::inserter( result, result.end() )
			) ;
			return result ;
		}
	}
}

#endif
