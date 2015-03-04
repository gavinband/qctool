
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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
#include "components/RelatednessComponent/mean_centre_genotypes.hpp"

// #define DEBUG_PCA_LOADING_COMPUTER 1

namespace {
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
			return string1 + genfile::string_utils::to_string( i+1 ) ;
		}
		else {
			return string2 + genfile::string_utils::to_string( i+1 -N ) ;
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
	m_sqrt_D_inverse = 1 / m_D.array().sqrt() ;
	m_U = udut.block( 0, 1, udut.rows(), n ) ;
	m_number_of_snps = number_of_snps ;
}

void PCALoadingComputer::begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) {
	assert( number_of_samples = std::size_t( m_U.rows() )) ;
	m_genotype_calls.resize( number_of_samples ) ;
	m_non_missingness.resize( number_of_samples ) ;
}

void PCALoadingComputer::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	data_reader.get(
		":genotypes:",
		genfile::vcf::get_threshholded_calls( m_genotype_calls, m_non_missingness, 0.9, 0, 0, 1, 2 )
	) ;
	assert( m_genotype_calls.size() == m_U.rows() ) ;
	assert( m_non_missingness.size() == m_U.rows() ) ;
	// setup the storage
	m_loading_vectors.resize( 2 * m_D.rows() ) ;
	m_loading_vectors.setConstant( std::numeric_limits< double >::quiet_NaN() ) ;
	double const allele_frequency = m_genotype_calls.sum() / ( 2.0 * m_non_missingness.sum() ) ;
	if( m_non_missingness.sum() > 0 && allele_frequency > 0.001 ) {
		//std::cerr << "pre-mean genotypes are: " << m_genotype_calls.head( 20 ).transpose() << "...\n" ;
		pca::mean_centre_genotypes( &m_genotype_calls, m_non_missingness, allele_frequency ) ;
		m_genotype_calls /= std::sqrt( 2.0 * allele_frequency * ( 1.0 - allele_frequency ) ) ;

#if DEBUG_PCA_LOADING_COMPUTER
		//std::cerr << "                    SNP: " << snp << ", allele frequency = " << allele_frequency << ".\n" ;
		//std::cerr << std::resetiosflags( std::ios::floatfield ) << std::setprecision( 5 ) ;
		//std::cerr << "pre-scale genotypes are: " << m_genotype_calls.head( 20 ).transpose() << "...\n" ;
		std::cerr << "          genotypes are: " << m_genotype_calls.head( 20 ).transpose() << "...\n" ;
		std::cerr << "     non-missingness is: " << m_non_missingness.head( 20 ).transpose() << "...\n" ;
		std::cerr << "                   U is: " << m_U.block(0,0,10,10) << "...\n" ;
		std::cerr << "                   D is: " << m_D << "...\n" ;
#endif // DEBUG_PCA_LOADING_COMPUTER
	

		//
		// Let X  be the L\times n matrix (L SNPs, n samples) of (mean-centred, scaled) genotypes.  We want
		// to compute the row of the matrix S of unit eigenvectors of the variance-covariance matrix
		// (1/L) X X^t that corresponds to the current SNP.
		// The matrix S is given by
		//               
		//       S = (1/√L) X U D^{-½}
		//
		// where
		//       (1/L) X^t X = U D U^t
		// is the eigenvalue decomposition of (1/L) X^t X that we are passed in via set_UDUT (and L is the number of SNPs).
		//
		// This is true since then
		//
		// S^t S = D^{-½} U^t (1/L) X^t X U D^{-½} = id
		//
		// (so columns of S are orthogonal) while
		//
		// (1/L X X^t) S = (1/L√L) X X^t X U D^{-½}
		//           = (1/√L) X U D U^t U D^{-½}
		//           = (1/√L) X U D^½
		//           = SD
		//
		// (so columns of S are eigenvectors with eigenvalues given by D.)
		//
#if 0
		m_loading_vectors.segment( 0, m_U.cols() ) =
			( m_genotype_calls.transpose() * m_U ) * m_D.array().sqrt().matrix().asDiagonal()
			/ ( ( m_D.transpose().array() * m_number_of_snps ).sqrt() ) ;
#else
		m_loading_vectors.segment( 0, m_U.cols() ) =
			( m_genotype_calls.transpose() * m_U ) * m_sqrt_D_inverse.asDiagonal() ;
		m_loading_vectors /= std::sqrt( m_number_of_snps ) ;
#endif
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
	}
	send_results(
		snp,
		m_non_missingness.sum(),
		allele_frequency,
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

void PCALoadingComputer::send_results( genfile::SNPIdentifyingData const& snp, double const N, double const frequency, Eigen::VectorXd const& data, GetNames get_names ) {
	m_result_signal( snp, N, frequency, data, get_names ) ;
}

std::string PCALoadingComputer::get_metadata() const {
	using namespace genfile::string_utils ;
	return "Number of SNPs: " + to_string( m_number_of_snps ) + "\n"
		+ "Number of samples: " + to_string( m_U.rows() ) + "\n"
		+ "These loadings represent unit eigenvectors of the variance-covariance matrix\n"
		+ "    1/L X X^t\n"
		+ "where X is the LxN matrix of genotypes at L SNPs and N samples (normalised across rows.)" ;
}

void PCALoadingComputer::end_processing_snps() {}



