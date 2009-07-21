#ifndef RANDOM_PERMUTATION_HPP
#define RANDOM_PERMUTATION_HPP

#include <iterator>
#include <algorithm>

// Randomly permute the given range.
// The random number generator provided must support the syntax
// random_number_generator( X )
// returning a number in the range [0, X] (where X >0 is an integer.)
template <typename Iterator, typename RNG>
void randomly_permute_range( Iterator begin, Iterator end, RNG const& random_number_generator ) {
	std::size_t N = std::distance( begin, end ) ;
	for( std::size_t i = 1; i <= N; ++i ) {
		std::size_t k = random_number_generator( N - i );
		if( k != ( N - i )) {
			std::swap( *(begin + k), *(begin + N - i)) ;
		}
	}
	
}

#endif
