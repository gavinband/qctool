#ifndef QCTOOL_RELATOTRON_HPP
#define QCTOOL_RELATOTRON_HPP

#include <boost/numeric/ublas/matrix.hpp>

#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/CohortIndividualSource.hpp"

#include "worker/Worker.hpp"

#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"

struct Relatotron: public genfile::SNPDataSourceProcessor::Callback
{
	static void declare_options( appcontext::OptionProcessor& options ) ;

	typedef genfile::SingleSNPGenotypeProbabilities SingleSNPGenotypeProbabilities ;
	typedef genfile::SNPIdentifyingData SNPIdentifyingData ;

	Relatotron( appcontext::OptionProcessor const& options, genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context ) ;
	
	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) ;
	void processed_snp( SNPIdentifyingData const& id_data, SingleSNPGenotypeProbabilities const& genotypes ) ;
	void end_processing_snps() ;
	
	// Run the relatedness algorithm
	void process( worker::Worker* worker = 0 ) ;
	
private:
	
	appcontext::OptionProcessor const& m_options ;
	genfile::CohortIndividualSource const& m_samples ;
	appcontext::UIContext& m_ui_context ;
	std::size_t m_number_of_samples ;
	double const m_probability_of_genotyping_error_per_snp ;
	
	std::vector< SNPIdentifyingData > m_snps ;
	std::vector< SingleSNPGenotypeProbabilities > m_genotypes ;
	std::vector< double > m_allele_frequencies ;
	
	typedef boost::numeric::ublas::vector< double > Vector ;
	typedef boost::numeric::ublas::matrix< double > Matrix ;
	typedef boost::numeric::ublas::zero_matrix< double > ZeroMatrix ;
	typedef boost::numeric::ublas::scalar_matrix< double > ConstantMatrix ;

	enum GenotypePair { e00 = 0, e01 = 1, e02 = 2, e10 = 3, e11 = 4, e12 = 5, e20 = 6, e21 = 7, e22 = 8 } ;
	
	std::vector< Matrix > m_genotype_per_ibd_matrices ;
	
	std::vector< std::size_t > m_row_samples ;
	std::vector< std::size_t > m_column_samples ;
	
private:
	
	void process_singlethreaded(
		Vector const& null_model,
		Vector const& alternative_model,
		Matrix* bf_matrix,
		std::vector< std::size_t > const& row_samples,
		std::vector< std::size_t > const& column_samples
	) ;

	void process_multithreaded(
		Vector const& null_model,
		Vector const& alternative_model,
		Matrix* bf_matrix,
		std::vector< std::size_t > const& row_samples,
		std::vector< std::size_t > const& column_samples,
		worker::Worker& worker
	) ;
	
	std::vector< std::size_t > parse_row_spec( std::string const& spec ) const ;
	
	Vector get_model_probabilities( std::vector< double > const& probabilities ) const ;
	
	void print_matrix( Matrix const& matrix, std::size_t const max_rows = 100, std::size_t const max_cols = 10 ) const ;
	double compute_maximum_likelihood_allele_frequency( SingleSNPGenotypeProbabilities const& genotypes ) const ;
	Matrix compute_genotype_probability_matrix( double const allele_frequency ) const ;

	void compute_pairwise_relatedness_log_bayes_factors(
		Vector const& null_ibd_probabilities,
		Vector const& alternative_ibd_probabilities,
		Matrix* result,
		std::vector< std::size_t > const& sample1_choice,
		std::vector< std::size_t > const& sample2_choice,
		appcontext::UIContext::ProgressContext* progress_context
	) const ;
	
	double compute_pairwise_relatedness_log_probability(
		std::size_t const sample1,
		std::size_t const sample2,
		Vector const& ibd_state_probabilities
	) const ;
	
	double compute_pairwise_relatedness_coefficients(
		std::size_t snp_i,
		std::size_t sample1,
		std::size_t sample2,
		Vector const& ibd_state_probabilities
	) const ;
	
	Matrix compute_probability_of_pair_of_genotypes_given_observations(
		std::size_t snp_i,
		std::size_t sample1,
		std::size_t sample2,
		double const theta // allele frequency
	) const ;
	
	double compute_probability_of_genotype_given_observations(
		std::size_t snp_i,
		std::size_t sample,
		std::size_t g,
		double const theta // allele frequency
	) const ;
	
	void write_relatedness_matrix(
		Matrix const& bf_matrix,
		std::string const& filename,
		std::vector< std::size_t > const& row_samples,
		std::vector< std::size_t > const& column_samples
	) const ;
} ;

#endif
