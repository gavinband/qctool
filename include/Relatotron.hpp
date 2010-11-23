#ifndef QCTOOL_RELATOTRON_HPP
#define QCTOOL_RELATOTRON_HPP


#include <limits>

#include <boost/numeric/ublas/matrix.hpp>

#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"

#include "fputils/floating_point_utils.hpp"

#include "appcontext/UIContext.hpp"

#include "OstreamTee.hpp"

struct Relatotron: public genfile::SNPDataSourceProcessor::Callback
{
	typedef genfile::SingleSNPGenotypeProbabilities SingleSNPGenotypeProbabilities ;
	typedef genfile::SNPIdentifyingData SNPIdentifyingData ;

	Relatotron( appcontext::UIContext& ui_context ):
	 	m_ui_context( ui_context )
	{} ;
	
	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
		m_number_of_samples = number_of_samples ;
		m_snps.clear() ;
		m_snps.reserve( number_of_snps ) ;
		m_genotypes.clear() ;
		m_genotypes.reserve( number_of_snps ) ;
		m_allele_frequencies.clear() ;
		m_allele_frequencies.reserve( number_of_snps ) ;
		m_genotype_per_ibd_matrices.clear() ;
		m_genotype_per_ibd_matrices.reserve( number_of_snps ) ;
	}

	void processed_snp( SNPIdentifyingData const& id_data, SingleSNPGenotypeProbabilities const& genotypes ) {
		assert( genotypes.get_number_of_samples() == m_number_of_samples ) ;
		assert( m_snps.size() == m_genotypes.size() ) ;
		m_snps.push_back( id_data ) ;
		m_genotypes.push_back( genotypes ) ;
		m_allele_frequencies.push_back( compute_maximum_likelihood_allele_frequency( genotypes )) ;
		m_genotype_per_ibd_matrices.push_back( compute_genotype_probability_matrix( m_allele_frequencies.back() )) ;
	}

	void end_processing_snps() {
		assert( m_snps.size() == m_genotypes.size() ) ;
		assert( m_snps.size() == m_allele_frequencies.size() ) ;

		std::size_t estimated_memory_usage = ( m_snps.capacity() * sizeof( SNPIdentifyingData ) )
			+ ( m_genotypes.capacity() * sizeof( SingleSNPGenotypeProbabilities ))
			+ ( m_genotypes.capacity() * m_number_of_samples * 3 * sizeof( double ))
			+ m_allele_frequencies.size() * sizeof( double )
			+ m_genotype_per_ibd_matrices.size() * sizeof( Matrix )
			+ m_genotype_per_ibd_matrices.size() * 27 * sizeof( double ) ;
		
		m_ui_context.logger()
			<< "Relatotron: finished loading "
			<< m_snps.size()
			<< " SNPs, memory usage is "
			<< std::fixed << std::setprecision( 1 ) << ( estimated_memory_usage / 1000000.0 )
			<< "Mb.\n" ;
	}
	
	void process() {
		Vector unrelated( 3 ) ;
		unrelated(0) = 1.0 ;
		unrelated(1) = 0.0 ;
		unrelated(2) = 0.0 ;

		Vector duplicate( 3 ) ;
		unrelated(0) = 0.0 ;
		unrelated(1) = 0.0 ;
		unrelated(2) = 1.0 ;
		
		Matrix bf_matrix( m_number_of_samples, m_number_of_samples ) ; // inefficient; this could be upper-triangular.
		{
			appcontext::UIContext::ProgressContext progress_context = m_ui_context.get_progress_context( "Calculating relatedness Bayes factors" ) ;
			compute_pairwise_relatedness_log_bayes_factors( unrelated, duplicate, &bf_matrix, progress_context ) ;
		}
		
		m_ui_context.logger() << "Top left of relatedness matrix: [\n" ;
		for( std::size_t i = 0; i < 10; ++i ) {
			for( std::size_t j = 0; j < 10; ++j ) {
				m_ui_context.logger() << std::setprecision( 2 ) << std::setw(5) << bf_matrix(i,j)  << " " ;
			}
			m_ui_context.logger()  << "\n" ;
		}
		m_ui_context.logger() << "]\n" ;
	}
	
private:
	
	appcontext::UIContext& m_ui_context ;
	std::size_t m_number_of_samples ;
	
	std::vector< SNPIdentifyingData > m_snps ;
	std::vector< SingleSNPGenotypeProbabilities > m_genotypes ;
	std::vector< double > m_allele_frequencies ;
	
	typedef boost::numeric::ublas::vector< double > Vector ;
	typedef boost::numeric::ublas::matrix< double > Matrix ;

	enum GenotypePair { e00 = 0, e01 = 1, e02 = 2, e10 = 3, e11 = 4, e12 = 5, e20 = 6, e21 = 7, e22 = 8 } ;
	
	std::vector< Matrix > m_genotype_per_ibd_matrices ;
	
private:
	double compute_maximum_likelihood_allele_frequency( SingleSNPGenotypeProbabilities const& genotypes ) const {
		// Under a model in which alleles are drawn randomly from a population of haplotypes,
		// with given allele frequency, return the maximum likelihood allele frequency.
		// Note that this does deal with null genotype calls (= missing data).
		double allele_count = 0.0 ;
		double data_count = 0.0 ;
		for( std::size_t i = 0; i < genotypes.size(); ++i ) {
			allele_count += genotypes.AB( i ) + 2.0 * genotypes.BB( i ) ; 
			data_count += genotypes.sum( i ) ;
		}

		if( data_count > 0.0 ) {
			return allele_count / 2 * data_count ;
		}
		else {
			return std::numeric_limits< double >::quiet_NaN() ;
		}
	}

	Matrix compute_genotype_probability_matrix(
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

	void compute_pairwise_relatedness_log_bayes_factors(
		Vector const& null_ibd_probabilities,
		Vector const& alternative_ibd_probabilities,
		Matrix* result,
		appcontext::UIContext::ProgressContext& progress_context
	) const {
		assert( result->size1() == m_number_of_samples ) ;
		assert( result->size2() == m_number_of_samples ) ;
		progress_context.notify_progress( 0, m_number_of_samples ) ;
		for( std::size_t sample1 = 0; sample1 < m_number_of_samples; ++sample1 ) {
			for( std::size_t sample2 = sample1; sample2 < m_number_of_samples; ++sample2 ) {
				(*result)( sample1, sample2 ) = compute_pairwise_relatedness_log_probability(
					sample1,
					sample2,
					alternative_ibd_probabilities
				)
				- compute_pairwise_relatedness_log_probability(
					sample1,
					sample2,
					null_ibd_probabilities
				) ;
			}
			progress_context.notify_progress( sample1 + 1, m_number_of_samples ) ;
		}
	}
	
	double compute_pairwise_relatedness_log_probability(
		std::size_t const sample1,
		std::size_t const sample2,
		Vector const& ibd_state_probabilities
	) const {
		std::vector< double > per_snp_log_probabilities( m_snps.size() ) ;
		for( std::size_t snp_i = 0; snp_i < m_snps.size(); ++snp_i ) {
			per_snp_log_probabilities[ snp_i ] = std::log(
				compute_pairwise_relatedness_coefficients(
					snp_i,
					sample1,
					sample2,
					ibd_state_probabilities
				)
			) ;
		}
		return fputils::log_sum_exp( per_snp_log_probabilities ) ;
	}
	
	double compute_pairwise_relatedness_coefficients(
		std::size_t snp_i,
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
					snp_i,
					sample1,
					sample2,
					m_allele_frequencies[ snp_i ]
				),
				m_genotype_per_ibd_matrices[ snp_i ]
			),
			ibd_state_probabilities
		)( 0 ) ;
	}
	
	Matrix compute_probability_of_pair_of_genotypes_given_observations(
		std::size_t snp_i,
		std::size_t sample1,
		std::size_t sample2,
		double const theta // allele frequency
	) const {
		Matrix result( 1, 9 ) ;
		for( std::size_t g1 = 0; g1 < 2; ++g1 ) {
			for( std::size_t g2 = 0; g2 < 2; ++g2 ) {
				result( 0, (3*g1)+g2 ) = compute_probability_of_genotype_given_observations(
					snp_i,
					sample1,
					g1,
					theta
				) *
				compute_probability_of_genotype_given_observations(
					snp_i,
					sample2,
					g2,
					theta
				) ;
			}
		}

		return result ;
	}

	double compute_probability_of_genotype_given_observations(
		std::size_t snp_i,
		std::size_t sample,
		std::size_t g,
		double const theta // allele frequency
	) const {
		double result = m_genotypes[ snp_i ]( sample, g ) ;
		double const null_call = m_genotypes[ snp_i ].null_call( sample ) ;
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
	
} ;

#endif
