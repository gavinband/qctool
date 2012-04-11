
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef RELATEDNESS_COMPONENT_MEAN_CENTRE_GENOTYPES_HPP
#define RELATEDNESS_COMPONENT_MEAN_CENTRE_GENOTYPES_HPP

#include <Eigen/Core>

namespace pca {
	void mean_centre_genotypes( 
		Eigen::VectorXd* threshholded_genotypes,
		Eigen::VectorXd const& non_missingness_matrix,
		double allele_frequency
	) ;
}

#endif
