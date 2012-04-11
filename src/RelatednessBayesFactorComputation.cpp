
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <limits>
#include <numeric>
#include <string>
#include <sstream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/bind.hpp>

#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/CohortIndividualSource.hpp"

#include "fputils/floating_point_utils.hpp"

#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "appcontext/OstreamTee.hpp"

#include "worker/Worker.hpp"
#include "worker/FunctionTask.hpp"

#include "string_utils/parse_utils.hpp"

#include "SampleBySampleComputation.hpp"
#include "RelatednessBayesFactorComputation.hpp"
#include "compute_maximum_likelihood_allele_frequency.hpp"

namespace impl {
	void print_matrix(
		appcontext::UIContext& ui_context,
		boost::numeric::ublas::matrix< double > const& matrix,
		std::size_t const max_rows = 10 ,
		std::size_t const max_cols = 10
	) {
		ui_context.logger() << "[\n" ;
		for( std::size_t i = 0; i < std::min( max_rows, matrix.size1() ); ++i ) {
			for( std::size_t j = 0; j < std::min( max_cols, matrix.size2() ); ++j ) {
				ui_context.logger() << std::setprecision( 3 ) << std::setw(8) << matrix( i, j ) << " " ;
			}
			ui_context.logger() << "\n" ;
		}
		ui_context.logger() << "].\n" ;
	}
}
RelatednessBayesFactorComputation::RelatednessBayesFactorComputation(
	appcontext::OptionProcessor const& options,
	appcontext::UIContext& ui_context
):
 	m_options( options ),
	m_ui_context( ui_context ),
	m_probability_of_genotyping_error_per_snp( options.get_value< double >( "-relatedness-epsilon" )),
	m_null_model_probabilities( get_model_probabilities( options.get_values< double >( "-relatedness-null" ))),
	m_alternative_model_probabilities( get_model_probabilities( options.get_values< double >( "-relatedness-alternative" )))
{}

std::string RelatednessBayesFactorComputation::get_summary() const {
	std::ostringstream ostr ;
	ostr << "Bayes factor, a la Thompson (1975) comparing models of relatedness: "
		<< m_null_model_probabilities << " (null) and "
		<< m_alternative_model_probabilities << " (alternative)." ;
	return ostr.str() ;
}

void RelatednessBayesFactorComputation::prepare(
	std::vector< SNPIdentifyingData > const& snps,
	std::vector< SingleSNPGenotypeProbabilities > const& genotypes
) {
	assert( snps.size() == genotypes.size() ) ;
	m_allele_frequencies.clear() ;
	m_allele_frequencies.reserve( snps.size() ) ;
	m_genotype_per_ibd_matrices.clear() ;
	m_genotype_per_ibd_matrices.reserve( snps.size() ) ;

	for( std::size_t snp_i = 0; snp_i < snps.size(); ++snp_i ) {
		m_allele_frequencies.push_back( compute_maximum_likelihood_allele_frequency( genotypes[ snp_i ] )) ;
		if( m_allele_frequencies.back() == m_allele_frequencies.back() ) { // if is not NaN 
			assert( m_allele_frequencies.back() >= 0.0 ) ;
			assert( m_allele_frequencies.back() <= 1.0 ) ;
			m_genotype_per_ibd_matrices.push_back( compute_genotype_probability_matrix( m_allele_frequencies.back() )) ;
		}
		else {
			m_genotype_per_ibd_matrices.push_back( Matrix( 0, 0 )) ;
			m_ui_context.logger() << "Allele frequency at SNP " << m_genotype_per_ibd_matrices.size() - 1 << " was not estimated.\n" ;
		}
	}
	
	assert( m_genotype_per_ibd_matrices.size() == m_allele_frequencies.size() ) ;
	m_ui_context.logger()
		<< "RelatednessBayesFactorComputation: first few allele frequencies are:\n" ;
	for( std::size_t i = 0; i < std::min( std::size_t( 10 ), m_allele_frequencies.size() ); ++i ) {
		m_ui_context.logger() << std::setprecision( 5 ) << std::setw( 8 ) << m_allele_frequencies[i] << " " ;
	}
	m_ui_context.logger() << "\n" ;

	m_ui_context.logger() << "RelatednessBayesFactorComputation: first few genotype-by-IBD matrices are:\n" ;
	for( std::size_t i = 0; i < std::min( m_genotype_per_ibd_matrices.size(), std::size_t( 3 )); ++i ) {
		m_ui_context.logger() << "SNP " << i << " (frequency = " << m_allele_frequencies[i] << "): " ;
		impl::print_matrix( m_ui_context, m_genotype_per_ibd_matrices[i] ) ;
	}
	m_ui_context.logger() << "\n" ;
}

RelatednessBayesFactorComputation::Vector RelatednessBayesFactorComputation::get_model_probabilities( std::vector< double > const& probabilities ) const {
	assert( probabilities.size() == 3 ) ;
	Vector result( 3 ) ;
	double sum = 0.0 ;
	for( std::size_t i = 0; i < 3; ++i ) {
		result(i) = probabilities[i] ;
		sum += probabilities[i] ;
	}
	assert( sum == 1.0 ) ;
	return result ;
}

RelatednessBayesFactorComputation::Matrix RelatednessBayesFactorComputation::compute_genotype_probability_matrix(
	double const allele_frequency
) const {
	// This matrix has 9 rows corresponding to the 3*3 possible genotypes of two samples,
	// and 3 columns corresponding to IBD0, IBD1, IBD2.
	Matrix result( 9, 3 ) ;

	double const pb = allele_frequency ; // shorthand
	double const pb_2 = std::pow( pb, 2.0 ) ;
	double const pb_3 = std::pow( pb, 3.0 ) ;
	double const pb_4 = std::pow( pb, 4.0 ) ;
	double const pa = 1.0 - allele_frequency ;
	double const pa_2 = std::pow( pa, 2.0 ) ;
	double const pa_3 = std::pow( pa, 3.0 ) ;
	double const pa_4 = std::pow( pa, 4.0 ) ;

	using std::pow ;

	// 0, 0
	result( e00, 0 ) = pa_4 ;
	result( e00, 1 ) = pa_3 ;
	result( e00, 2 ) = pa_2 ;
	// 0, 1
	result( e01, 0 ) = 2.0*pa_3*pb ;
	result( e01, 1 ) = pa_2*pb ;
	result( e01, 2 ) = 0.0 ;
	// 0, 2
	result( e02, 0 ) = pa_2 * pb_2 ;
	result( e02, 1 ) = 0.0 ;
	result( e02, 2 ) = 0.0 ;
	// 1, 0
	result( e10, 0 ) = 2.0*pa_3*pb ;
	result( e10, 1 ) = pa_2 * pb ;
	result( e10, 2 ) = 0.0 ;
	// 1, 1
	result( e11, 0 ) = 4.0*pa_2*pb_2 ;
	result( e11, 1 ) = pa*pb*(pa+pb) ;
	result( e11, 2 ) = 2.0*pa*pb ;
	// 1, 2
	result( e12, 0 ) = 2*pa*pb_3 ;
	result( e12, 1 ) = pa*pb_2 ;
	result( e12, 2 ) = 0.0 ;
	// 2, 0
	result( e20, 0 ) = pa_2*pb_2 ;
	result( e20, 1 ) = 0.0 ;
	result( e20, 2 ) = 0.0 ;
	// 2, 1
	result( e21, 0 ) = 2*pa*pb_3 ;
	result( e21, 1 ) = pa*pb_2 ;
	result( e21, 2 ) = 0.0 ;
	// 2, 2
	result( e22, 0 ) = pb_4 ;
	result( e22, 1 ) = pb_3 ;
	result( e22, 2 ) = pb_2 ;

	return result ;
}

double RelatednessBayesFactorComputation::operator()(
	std::size_t const sample1,
	std::size_t const sample2,
	std::vector< SingleSNPGenotypeProbabilities > const& genotypes
) {
	return compute_pairwise_relatedness_log_probability( sample1, sample2, genotypes, m_alternative_model_probabilities )
		- compute_pairwise_relatedness_log_probability( sample1, sample2, genotypes, m_null_model_probabilities ) ;
}

double RelatednessBayesFactorComputation::compute_pairwise_relatedness_log_probability(
	std::size_t const sample1,
	std::size_t const sample2,
	std::vector< SingleSNPGenotypeProbabilities > const& genotypes,
	Vector const& ibd_state_probabilities
) const {
	std::vector< double > per_snp_log_probabilities( genotypes.size() ) ;
	for( std::size_t snp_i = 0; snp_i < genotypes.size(); ++snp_i ) {
		if( m_allele_frequencies[ snp_i ] == m_allele_frequencies[ snp_i ] ) { // if not NaN
			// To be tolerant to genotyping error, with probability m_genotype_error_probability,
			// we ignore the SNP's data.  Otherwise we use it.
			per_snp_log_probabilities[ snp_i ] = std::log(
				(1 - m_probability_of_genotyping_error_per_snp ) * compute_pairwise_relatedness_coefficients(
					genotypes[ snp_i ],
					m_allele_frequencies[ snp_i ],
					m_genotype_per_ibd_matrices[ snp_i ],
					sample1,
					sample2,
					ibd_state_probabilities
				)
				+
				m_probability_of_genotyping_error_per_snp
			) ;
			//std::cerr << "probability ( " << ibd_state_probabilities << " ) for " << sample1 << ", " << sample2 << ", snp " << snp_i << " is " << std::exp( per_snp_log_probabilities[ snp_i ] ) << ".\n" ;
		}
		else {
			// ignore the snp.
			per_snp_log_probabilities[ snp_i ] = 0.0 ;
		}
	}
	double result = std::accumulate( per_snp_log_probabilities.begin(), per_snp_log_probabilities.end(), 0.0 ) ;
	//std::cerr << "model " << ibd_state_probabilities << ": likelihood for " << sample1 << ", " << sample2 << " is " << std::exp( result ) << ".\n" ;
	return result ;
}

double RelatednessBayesFactorComputation::compute_pairwise_relatedness_coefficients(
	SingleSNPGenotypeProbabilities const& genotype,
	double allele_frequency,
	Matrix const& genotype_per_ibd_matrix,
	std::size_t sample1,
	std::size_t sample2,
	Vector const& ibd_state_probabilities
) const {
	// Compute the product of the vector of genotype-pair probabilities
	// times the genotype-per-IBD matrix.
	// Then 
	using namespace boost::numeric::ublas ;
	return prod(
		prod(
			compute_probability_of_pair_of_genotypes_given_observations(
				genotype,
				sample1,
				sample2,
				allele_frequency
			),
			genotype_per_ibd_matrix
		),
		ibd_state_probabilities
	)( 0 ) ;
}

RelatednessBayesFactorComputation::Matrix RelatednessBayesFactorComputation::compute_probability_of_pair_of_genotypes_given_observations(
	SingleSNPGenotypeProbabilities const& genotype,
	std::size_t sample1,
	std::size_t sample2,
	double const theta // allele frequency
) const {
	Matrix result( 1, 9 ) ;

	for( std::size_t g1 = 0; g1 < 3; ++g1 ) {
		double p1 = compute_probability_of_genotype_given_observations(
			genotype,
			sample1,
			g1,
			theta
		) ;
		result( 0, (3*g1) ) = p1 ;
		result( 0, (3*g1) + 1 ) = p1 ;
		result( 0, (3*g1) + 2 ) = p1 ;
	}
	
	for( std::size_t g2 = 0; g2 < 3; ++g2 ) {
		double p2 = compute_probability_of_genotype_given_observations(
			genotype,
			sample2,
			g2,
			theta
		) ;
		result( 0, g2 ) *= p2 ;
		result( 0, 3+g2 ) *= p2 ;
		result( 0, 6+g2 ) *= p2 ;
	}
	/*
	std::cerr << "compute_probability_of_pair_of_genotypes_given_observations( "
		<< snp_i << ", " << sample1 << ", " << sample2 << ", " << theta << " ): result is:\n" ;
	std::cerr << result << "\n" ;
	*/	
	return result ;
}

double RelatednessBayesFactorComputation::compute_probability_of_genotype_given_observations(
	SingleSNPGenotypeProbabilities const& genotype,
	std::size_t sample,
	std::size_t g,
	double const theta // allele frequency
) const {
	double result = genotype( sample, g ) ;
	double const null_call = genotype.null_call( sample ) ;
	// If  genotype is missing (null call), we fill in the probability
	// from the allele frequency.
	switch( g ) {
		case 0:
			result += null_call * (1-theta) * (1-theta) ;
			break ;
		case 1:
			result += null_call * 2.0 * theta * ( 1 - theta ) ;
			break ;
		case 2:
			result += null_call * theta * theta ;
			break ;
		default:
			assert(0) ;
	}
	
	return result ;
}
