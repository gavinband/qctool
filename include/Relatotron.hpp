#ifndef QCTOOL_RELATOTRON_HPP
#define QCTOOL_RELATOTRON_HPP


#include <limits>
#include <numeric>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/CohortIndividualSource.hpp"

#include "fputils/floating_point_utils.hpp"

#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/get_current_time_as_string.hpp"

#include "string_utils/parse_utils.hpp"

#include "OstreamTee.hpp"

struct Relatotron: public genfile::SNPDataSourceProcessor::Callback
{
	static void declare_options( appcontext::OptionProcessor& options ) {
		options.declare_group( "Relatedness options" ) ;
		options[ "-relatedness" ]
			.set_description( "Compute relatedness matrices pairwise for all samples.  (This can take a long time)." )
			.set_takes_single_value() ;
		options[ "-relatedness-samples" ]
			.set_description( "Choose ranges of samples to compute relatedness for."
				" This option should be in the form a-b,c-d where a-b and c-d are the sample ranges between which"
				" to compute relatedness." )
			.set_takes_single_value() ;
		options[ "-relatedness-epsilon" ]
			.set_description( "Set the probability of genotyping error at a SNP."
				" This is used to make the model tolerant to genotyping errors." )
			.set_takes_single_value()
			.set_default_value( 1/1000.0 ) ;
	}
	
	typedef genfile::SingleSNPGenotypeProbabilities SingleSNPGenotypeProbabilities ;
	typedef genfile::SNPIdentifyingData SNPIdentifyingData ;

	Relatotron( appcontext::OptionProcessor const& options, genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context ):
		m_options( options ),
		m_samples( samples ),
	 	m_ui_context( ui_context ),
		m_probability_of_genotyping_error_per_snp( options.get_value< double >( "-relatedness-epsilon" ))
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
		assert( m_allele_frequencies.back() >= 0.0 ) ;
		assert( m_allele_frequencies.back() <= 1.0 ) ;
		m_genotype_per_ibd_matrices.push_back( compute_genotype_probability_matrix( m_allele_frequencies.back() )) ;
		m_ui_context.logger() << "Genotype per IBD matrix (SNP " << m_genotype_per_ibd_matrices.size() - 1 << " is: [\n" ;
		print_matrix( m_genotype_per_ibd_matrices.back() ) ;
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
			
		m_ui_context.logger()
			<< "Relatotron: first few allele frequencies are:\n" ;
		for( std::size_t i = 0; i < std::min( std::size_t( 10 ), m_allele_frequencies.size() ); ++i ) {
			m_ui_context.logger() << std::setprecision( 5 ) << std::setw( 8 ) << m_allele_frequencies[i] << " " ;
		}
		m_ui_context.logger() << "\n" ;
	}
	
	void process() {
		Vector unrelated( 3 ) ;
		unrelated(0) = 1.0 ;
		unrelated(1) = 0.0 ;
		unrelated(2) = 0.0 ;

		Vector duplicate( 3 ) ;
		duplicate(0) = 0.0 ;
		duplicate(1) = 0.0 ;
		duplicate(2) = 1.0 ;
		
		// Store the results in a plain matrix.  This is inefficient storage-wise.
		Matrix bf_matrix = ConstantMatrix( m_number_of_samples, m_number_of_samples, -std::numeric_limits< double >::infinity() ) ; 
		
		if( m_options.check_if_option_was_supplied( "-relatedness-samples" )) {
			std::string const range_spec = m_options.get_value< std::string >( "-relatedness-samples" ) ;
			std::vector< std::string > range_specs = string_utils::split_and_strip( range_spec, "," ) ;
			if( !range_specs.size() == 2 ) {
				throw genfile::BadArgumentError( "Relatotron::process()", "value \"" + range_spec + "\" for option -relatedness-samples") ;
			}
			std::vector< std::pair< std::size_t, std::size_t > > ranges ;
			ranges.push_back( string_utils::parse_range< std::size_t >( range_specs[0] )) ;
			ranges.push_back( string_utils::parse_range< std::size_t >( range_specs[1] )) ;
			
			if( ranges.front().first > m_number_of_samples || ranges.front().second > m_number_of_samples ) {
				throw genfile::BadArgumentError( "Relatotron::process()", "first range spec exceeds number of samples (" + string_utils::to_string( m_number_of_samples ) + ")") ;
			}
			if( ranges.back().first > m_number_of_samples || ranges.back().second > m_number_of_samples ) {
				throw genfile::BadArgumentError( "Relatotron::process()", "second range spec exceeds number of samples (" + string_utils::to_string( m_number_of_samples ) + ")") ;
			}
			
			appcontext::UIContext::ProgressContext progress_context = m_ui_context.get_progress_context( "Calculating relatedness Bayes factors" ) ;
			compute_pairwise_relatedness_log_bayes_factors(
				unrelated,
				duplicate,
				&bf_matrix,
				ranges.front(),
				ranges.back(),
				progress_context
			) ;
		}
		else {
			appcontext::UIContext::ProgressContext progress_context = m_ui_context.get_progress_context( "Calculating relatedness Bayes factors" ) ;
			compute_pairwise_relatedness_log_bayes_factors( unrelated, duplicate, &bf_matrix, progress_context ) ;
		}
		
		m_ui_context.logger() << "Top left of relatedness matrix: [\n" ;
		for( std::size_t i = 0; i < std::min( std::size_t( 10 ), m_number_of_samples ); ++i ) {
			for( std::size_t j = 0; j < std::min( std::size_t( 10 ), m_number_of_samples ); ++j ) {
				m_ui_context.logger() << std::setprecision( 2 ) << std::setw(8) << std::exp( bf_matrix(i,j) )  << " " ;
			}
			m_ui_context.logger()  << "\n" ;
		}
		m_ui_context.logger() << "]\n" ;
		
		if( m_options.check_if_option_was_supplied( "-relatedness" )) {
			write_relatedness_matrix( bf_matrix, m_options.get_value< std::string >( "-relatedness" ) ) ;
		}
	}
	
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
	
private:
	
	void print_matrix( Matrix const& matrix, std::size_t const max_rows = 100, std::size_t const max_cols = 10 ) const {
		m_ui_context.logger() << "[\n" ;
		for( std::size_t i = 0; i < std::min( max_rows, matrix.size1() ); ++i ) {
			for( std::size_t j = 0; j < std::min( max_cols, matrix.size2() ); ++j ) {
				m_ui_context.logger() << std::setprecision( 3 ) << std::setw(8) << matrix( i, j ) << " " ;
			}
			m_ui_context.logger() << "\n" ;
		}
		m_ui_context.logger() << "].\n" ;
	}

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
			return allele_count / (2.0 * data_count ) ;
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
		return compute_pairwise_relatedness_log_bayes_factors(
			null_ibd_probabilities,
			alternative_ibd_probabilities,
			result,
			std::make_pair( std::size_t( 0 ), m_number_of_samples - 1 ),
			std::make_pair( std::size_t( 0 ), m_number_of_samples - 1 ),
			progress_context
		) ;
	}

	void compute_pairwise_relatedness_log_bayes_factors(
		Vector const& null_ibd_probabilities,
		Vector const& alternative_ibd_probabilities,
		Matrix* result,
		std::pair< std::size_t, std::size_t > const& sample1_range,
		std::pair< std::size_t, std::size_t > const& sample2_range,
		appcontext::UIContext::ProgressContext& progress_context
	) const {
		assert( result->size1() == m_number_of_samples ) ;
		assert( result->size2() == m_number_of_samples ) ;
		assert( sample1_range.first < m_number_of_samples ) ;
		assert( sample1_range.second < m_number_of_samples ) ;
		assert( sample2_range.first < m_number_of_samples ) ;
		assert( sample2_range.second < m_number_of_samples ) ;
		std::cerr << sample1_range.first << ":" << sample1_range.second << ", " << sample2_range.first << ":" << sample2_range.second << ".\n" ;
		progress_context.notify_progress( 0, m_number_of_samples ) ;
		for( std::size_t sample1 = sample1_range.first; sample1 <= sample1_range.second; ++sample1 ) {
			for( std::size_t sample2 = std::max( sample2_range.first, sample1 ); sample2 <= sample2_range.second; ++sample2 ) {
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
			// To be tolerant to genotyping error, with probability m_genotype_error_probability,
			// we ignore the SNP's data.  Otherwise we use it.
			per_snp_log_probabilities[ snp_i ] = std::log(
				(1 - m_probability_of_genotyping_error_per_snp ) * compute_pairwise_relatedness_coefficients(
					snp_i,
					sample1,
					sample2,
					ibd_state_probabilities
				)
				+
				m_probability_of_genotyping_error_per_snp
			) ;
			//std::cerr << "probability ( " << ibd_state_probabilities << " ) for " << sample1 << ", " << sample2 << ", snp " << snp_i << " is " << std::exp( per_snp_log_probabilities[ snp_i ] ) << ".\n" ;

		}
		double result = std::accumulate( per_snp_log_probabilities.begin(), per_snp_log_probabilities.end(), 0.0 ) ;
		//std::cerr << "model " << ibd_state_probabilities << ": likelihood for " << sample1 << ", " << sample2 << " is " << std::exp( result ) << ".\n" ;
		return result ;
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

		for( std::size_t g1 = 0; g1 < 3; ++g1 ) {
			double p1 = compute_probability_of_genotype_given_observations(
				snp_i,
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
				snp_i,
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
	
	void write_relatedness_matrix( Matrix const& bf_matrix, std::string const& filename ) const {
		appcontext::OUTPUT_FILE_PTR file = open_file_for_output( filename ) ;
		(*file) << "# Written by Relatotron, " << appcontext::get_current_time_as_string() << "\n" ;
		if( m_genotypes.size() > 0 ) {
			assert( bf_matrix.size1() == m_genotypes[0].size() ) ;
			assert( bf_matrix.size2() == m_genotypes[0].size() ) ;
			assert( bf_matrix.size2() == m_samples.get_number_of_individuals() ) ;
			
			(*file) << "id," ;
			for( std::size_t i = 0; i < m_samples.get_number_of_individuals(); ++i ) {
				if( i > 0 ) {
					(*file) << "," ;
				}
				(*file) << m_samples.get_entry( i, "id_1" ).as< std::string >() ;
			}
			(*file) << "\n" ;

			for( std::size_t i = 0; i < bf_matrix.size1(); ++i ) {
				(*file) << m_samples.get_entry( i, "id_1" ).as< std::string >() << "," ;
				for( std::size_t j = 0; j < bf_matrix.size2(); ++j ) {
					if( j > 0 ) {
						(*file) << "," ;
					}
					(*file) << bf_matrix( i, j ) ;
				}
				(*file) << "\n" ;
			}
		}
	}
	
} ;

#endif
