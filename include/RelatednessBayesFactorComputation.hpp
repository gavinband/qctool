#ifndef QCTOOL_RELATEDNESSBAYESFACTORCOMPUTATION_HPP
#define QCTOOL_RELATEDNESSBAYESFACTORCOMPUTATION_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/CohortIndividualSource.hpp"

#include "worker/Worker.hpp"

#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "SampleBySampleComputation.hpp"

struct RelatednessBayesFactorComputation: public SampleBySampleComputation
{
	static void declare_options( appcontext::OptionProcessor& options ) ;

	typedef genfile::SingleSNPGenotypeProbabilities SingleSNPGenotypeProbabilities ;
	typedef genfile::SNPIdentifyingData SNPIdentifyingData ;

	RelatednessBayesFactorComputation( appcontext::OptionProcessor const& options, appcontext::UIContext& ui_context ) ;
	
	void prepare( std::vector< SNPIdentifyingData > const& snps, std::vector< SingleSNPGenotypeProbabilities > const& genotypes ) ;
	
	// Run the relatedness algorithm
	double operator()(
		std::size_t const sample1,
		std::size_t const sample2,
		std::vector< SingleSNPGenotypeProbabilities > const& genotypes
	) ;
	
	std::string get_summary() const ;
	
private:
	
	typedef boost::numeric::ublas::vector< double > Vector ;
	typedef boost::numeric::ublas::matrix< double > Matrix ;
	typedef boost::numeric::ublas::zero_matrix< double > ZeroMatrix ;
	typedef boost::numeric::ublas::scalar_matrix< double > ConstantMatrix ;
	enum GenotypePair { e00 = 0, e01 = 1, e02 = 2, e10 = 3, e11 = 4, e12 = 5, e20 = 6, e21 = 7, e22 = 8 } ;

	appcontext::OptionProcessor const& m_options ;
	appcontext::UIContext& m_ui_context ;
	std::size_t m_number_of_samples ;
	double const m_probability_of_genotyping_error_per_snp ;
	Vector m_null_model_probabilities ;
	Vector m_alternative_model_probabilities ;
	
	std::vector< double > m_allele_frequencies ;
	std::vector< Matrix > m_genotype_per_ibd_matrices ;
private:
	
	Vector get_model_probabilities( std::vector< double > const& probabilities ) const ;
	
	Matrix compute_genotype_probability_matrix( double const allele_frequency ) const ;

	double compute_pairwise_relatedness_log_probability(
		std::size_t const sample1,
		std::size_t const sample2,
		std::vector< SingleSNPGenotypeProbabilities > const& genotypes,
		Vector const& ibd_state_probabilities
	) const ;
	
	double compute_pairwise_relatedness_coefficients(
		SingleSNPGenotypeProbabilities const& genotypes,
		double allele_frequency,
		Matrix const& genotype_per_ibd_matrix,
		std::size_t sample1,
		std::size_t sample2,
		Vector const& ibd_state_probabilities
	) const ;
	
	Matrix compute_probability_of_pair_of_genotypes_given_observations(
		SingleSNPGenotypeProbabilities const& genotypes,
		std::size_t sample1,
		std::size_t sample2,
		double const theta // allele frequency
	) const ;
	
	double compute_probability_of_genotype_given_observations(
		SingleSNPGenotypeProbabilities const& genotypes,
		std::size_t sample,
		std::size_t g,
		double const theta // allele frequency
	) const ;
} ;

#endif
