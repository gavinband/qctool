#include <iostream>
#include <iomanip>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/function.hpp>
#include "../config.hpp"
#if HAVE_CBLAS
	#include "cblas.h"
#endif
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "components/RelatednessComponent/PCALoadingComputer.hpp"
#include "components/RelatednessComponent/LapackEigenDecomposition.hpp"

namespace {
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
	
	template< typename Vector1, typename Vector2, typename NonMissingVector >
	double compute_correlation( Vector1 const& v1, Vector2 const& v2, NonMissingVector const& non_missingness_indicator ) {
		assert( v1.size() == v2.size() ) ;
		double non_missingness = non_missingness_indicator.sum() ;
		double mean1 = 0.0 ;
		double mean2 = 0.0 ;
		for( int i = 0; i < v1.size(); ++i ) {
			if( non_missingness_indicator( i )) {
				mean1 += v1(i) / non_missingness ;
				mean2 += v2(i) / non_missingness ;
			}
		}

		double covariance = 0.0 ;
		double variance1 = 0.0 ;
		double variance2 = 0.0 ;
		for( int i = 0; i < v1.size(); ++i ) {
			if( non_missingness_indicator( i )) {
				covariance += ( v1(i) - mean1 ) * ( v2(i) - mean2 ) ;
				variance1 += ( v1(i) - mean1 ) * ( v1(i) - mean1 ) ;
				variance2 += ( v2(i) - mean2 ) * ( v2(i) - mean2 ) ;
			}
		}
		
		// We should divide the covariance by N-1 and also
		// divide each variance by the same quantity.
		// But this washes out in the ratio.
		
		return covariance / std::sqrt( variance1 * variance2 ) ;
	}
	
	std::string eigenvector_column_names( std::size_t N, std::string const& string1, std::string const& string2, std::size_t i ) {
		if( i < N ) {
			return string1 + genfile::string_utils::to_string( i ) ;
		}
		else {
			return string2 + genfile::string_utils::to_string( i - N ) ;
		}
	}
}

PCALoadingComputer::PCALoadingComputer( int number_of_loadings ):
	m_number_of_loadings( number_of_loadings ),
	m_number_of_snps( 1 )
{}

void PCALoadingComputer::set_UDUT( std::size_t number_of_snps, Matrix const& udut ) {
	assert( udut.cols() == udut.rows() + 1 ) ;
	int n = std::min( int( m_number_of_loadings ), int( udut.rows() ) ) ;
	m_D = udut.block( 0, 0, n, 1 ) ;
	m_U = udut.block( 0, 1, udut.rows(), n ) ;
	m_number_of_snps = number_of_snps ;
}

void PCALoadingComputer::begin_processing_snps( std::size_t number_of_samples ) {
	assert( number_of_samples = std::size_t( m_U.rows() )) ;
	m_genotype_calls.resize( number_of_samples ) ;
	m_non_missingness.resize( number_of_samples ) ;
}

void PCALoadingComputer::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	genfile::vcf::ThreshholdingGenotypeSetter< Eigen::VectorXd > setter( m_genotype_calls, m_non_missingness, 0.9, 0, 0, 1, 2 ) ;
	data_reader.get( "genotypes", setter ) ;
	assert( m_genotype_calls.size() == m_U.rows() ) ;
	assert( m_non_missingness.size() == m_U.rows() ) ;
	double const allele_frequency = m_genotype_calls.sum() / ( 2.0 * m_non_missingness.sum() ) ;
	//std::cerr << "pre-mean genotypes are: " << m_genotype_calls.head( 20 ).transpose() << "...\n" ;
	mean_centre_genotypes( &m_genotype_calls, m_non_missingness, allele_frequency ) ;

	//std::cerr << "                    SNP: " << snp << ", allele frequency = " << allele_frequency << ".\n" ;
	//std::cerr << std::resetiosflags( std::ios::floatfield ) << std::setprecision( 5 ) ;
	//std::cerr << "pre-scale genotypes are: " << m_genotype_calls.head( 20 ).transpose() << "...\n" ;
	m_genotype_calls /= std::sqrt( 2.0 * allele_frequency * ( 1.0 - allele_frequency ) ) ;
	//std::cerr << "          genotypes are: " << m_genotype_calls.head( 20 ).transpose() << "...\n" ;
	//std::cerr << "     non-missingness is: " << m_non_missingness.head( 20 ).transpose() << "...\n" ;
	//std::cerr << "                   U is: " << m_U.block(0,0,10,10) << "...\n" ;
	//std::cerr << "                   D is: " << m_D << "...\n" ;

	//
	// Let X  be the L\times n matrix (L SNPs, n samples) of (mean-centred, scaled) genotypes.  We want
	// to compute the row of the matrix S of unit eigenvectors of X X^t that corresponds to the current SNP.
	// The matrix S is given by
	//               1 
	//       S = ------- X U D^{-1/2}
	//           \sqrt(L)
	// where
	//       (1/L) X^t X = U D U^t
	// is the eigenvalue decomposition of (1/L) X^t X that we are passed in via set_UDUT (and L is the number of SNPs).
	//
	// Here we must compute a single row of S.
	m_loading_vectors.resize( 2 * m_D.rows() ) ;

	m_loading_vectors.setConstant( std::numeric_limits< double >::quiet_NaN() ) ;
	m_loading_vectors.segment( 0, m_U.cols() ) = ( m_genotype_calls.transpose() * m_U ).array() / ( ( m_D.transpose().array() * m_number_of_snps ).sqrt() ) ;
	// We also wish to compute the correlation between the SNP and the PCA component.
	// With S as above, the PCA components are the projections of columns of X onto columns of S.
	// If we want samples to correspond to columns, this is
	//   S^t X 
	// which can be re-written
	//   sqrt(L) U D^{1/2}
	// i.e. we may as well compute the correlation with columns of U.
	if( m_non_missingness.sum() > 10 ) {
		for( int i = 0; i < m_U.cols(); ++i ) {
			m_loading_vectors( m_U.cols() + i ) = compute_correlation( m_genotype_calls, m_U.col( i ), m_non_missingness ) ;
		}
	}
	
	send_results(
		snp,
		m_loading_vectors,
		boost::bind(
			&eigenvector_column_names,
			m_U.cols(),
			"eigenvector_",
			"correlation_",
			_1
		)
	) ;
}

void PCALoadingComputer::send_results_to( ResultCallback callback ) {
	m_result_signal.connect( callback ) ;
}

void PCALoadingComputer::send_results( genfile::SNPIdentifyingData const& snp, Eigen::VectorXd const& data, GetNames get_names ) {
	m_result_signal( snp, data, get_names ) ;
}

std::string PCALoadingComputer::get_metadata() const {
	using namespace genfile::string_utils ;
	return "Number of SNPs: " + to_string( m_number_of_snps ) + "\n"
		+ "Number of samples: " + to_string( m_genotype_calls.size() ) ;
}

void PCALoadingComputer::end_processing_snps() {}



