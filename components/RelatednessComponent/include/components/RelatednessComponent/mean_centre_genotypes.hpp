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
