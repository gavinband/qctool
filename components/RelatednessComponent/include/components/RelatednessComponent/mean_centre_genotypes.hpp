
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef RELATEDNESS_COMPONENT_MEAN_CENTRE_GENOTYPES_HPP
#define RELATEDNESS_COMPONENT_MEAN_CENTRE_GENOTYPES_HPP

#include <Eigen/Core>

namespace pca {
	template< typename Vector1, typename Vector2 >
	void mean_centre_genotypes( 
		Vector1* threshholded_genotypes,
		Vector2 const& non_missingness_matrix,
		double allele_frequency
	) {
		std::size_t const number_of_samples = threshholded_genotypes->size() ;
		assert( std::size_t( non_missingness_matrix.size() ) == number_of_samples ) ;
		for( std::size_t sample_i = 0; sample_i < number_of_samples; ++sample_i ) {
			if( non_missingness_matrix( sample_i ) ) {
				(*threshholded_genotypes)( sample_i ) -= 2.0 * allele_frequency ;
			}
			else {
				(*threshholded_genotypes)( sample_i ) = 0.0 ; // this sample does not contribute for this SNP.
			}
		}
	}
}

#endif
