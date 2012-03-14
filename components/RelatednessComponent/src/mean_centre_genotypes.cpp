#include <Eigen/Core>

namespace pca {
	void mean_centre_genotypes( 
		Eigen::VectorXd* threshholded_genotypes,
		Eigen::VectorXd const& non_missingness_matrix,
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
