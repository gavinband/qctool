#ifndef QCTOOL_HAPLOTYPE_FREQUENCY_COMPONENT_HPP
#define QCTOOL_HAPLOTYPE_FREQUENCY_COMPONENT_HPP

#include <string>
#include <boost/function.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/RandomAccessSNPDataSource.hpp"

// Loglikelihood of table of genotypes at 2 SNPs, given parameters specifying the haplotype frequencies.
// The parameters are \pi_01, \pi_10, and \pi_11.  Then \pi_00 is one minus the sum of the others.
// \pi_ab is the frequency of the haplotype with genotype a at SNP 1 and b at SNP 2.
// Now the probability of any table is obtained
// 
struct HaplotypeFrequencyLogLikelihood {
	HaplotypeFrequencyLogLikelihood( Eigen::Matrix const& genotype_table ) ;
	void evaluate_at( Vector const& parameters ) ;
	double get_value_of_function() {
		double result = 0.0 ;
		for( std::size_t g1 = 0; g1 < 3; ++g1 ) {
			for( std::size_t g2 = 0; g2 < 3; ++g2 ) {
				
			}
		}
	}
	Vector get_value_of_first_derivative() ;
	Matrix get_value_of_second_derivative() ;
} ;

struct HaplotypeFrequencyComponent: public genfile::SNPDataSourceProcessor::Callback {
public:
	typedef std::auto_ptr< CallComparerComponent > UniquePtr ;

	static void declare_options( appcontext::OptionProcessor& options ) {
		options.declare_group( "LD computation options" ) ;
		options[ "-compute-ld" ]
			.set_description( "Compute LD pairwise metrics between the main dataset and SNPs." )
			.set_takes_single_value() ;
		options[ "-ld-source-snps" ]
			.set_description( "Specify the name of a data source containing SNPs with which LD should be computed." )
			.set_takes_single_value() ;
		options.option_implies_option( "-compute-ld", "-ld-source-snps" ) ;
		options.option_implies_option( "-ld-source-snps", "-compute-ld" ) ;
	}
	
	static UniquePtr create( appcontext::OptionProcessor const& options ) {
		HaplotypeFrequencyComponent::UniquePtr result ;
		result.reset(
			new HaplotypeFrequencyComponent(
				genfile::SNPDataSourceChain::create(
					genfile::wildcard::find_files_by_chromosome(
						options[ "-compute-ld" ],
						genfile::wildcard::e_AUTOSOMAL_CHROMOSOMES
					)
				)
			)
		) ;
		
		
	}

public:
	HaplotypeFrequencyComponent(
		genfile::SNPDataSource::UniquePtr source
	):
		m_source( source ),
		m_threshhold( 0.9 )
	{}

public:
	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
		assert( m_source->number_of_samples() == number_of_samples () ) ;
	}

	void processed_snp( genfile::SNPIdentifyingData const& target_snp, genfile::VariantDataReader& target_data_reader ) {
		genfile::SNPIdentifyingData source_snp ;
		genfile::SingleSNPGenotypeProbabilities source_probs, target_probs ;
		data_reader.get( "genotypes", target_probs ) ;
		m_source->reset_to_start() ;
		while( m_source->get_snp_identifying_data( source_snp )) {
			genfile::VariantDataReader::UniquePtr source_data_reader = m_source->read_variant_data() ;
			compute_ld( source_snp, target_snp, *source_data_reader, target_data_reader ) ;
		}
	}
	
	void compute_ld_measures(
		genfile::SNPIdentifyingData const& source_snp,
		genfile::VariantDataReader& source_data_reader,
		genfile::SNPIdentifyingData const& target_snp,
		genfile::VariantDataReader& target_data_reader
	) {
		std::vector< std::vector< int > > genotypes( 2 ) ;
		GenotypeSetter< std::vector< int > > source_getter( genotypes[0], m_threshhold ) ;
		GenotypeSetter< std::vector< int > > target_getter( genotypes[1], m_threshhold ) ;
		source_data_reader.get( "genotypes", source_getter ) ;
		source_data_reader.get( "genotypes", target_getter ) ;
		assert( genotypes[0].size() == m_source->number_of_samples() ) ;
		assert( genotypes[0].size() == genotypes[1].size() ) ;
		compute_ld(
			source_snp,
			target_snp,
			genotypes
		) ;
	}
	
	void compute_ld_measures(
		genfile::SNPIdentifyingData const& source_snp,
		genfile::SNPIdentifyingData const& target_snp,
		std::vector< std::vector< int > > const& genotypes
	) {
		// Construct table of genotypes at each SNP.
		Eigen::Matrix3d table = Eigen::Matrix3d::Zero() ;
		for( std::size_t i = 0; i < m_source->number_of_samples(); ++i ) {
			if( genotypes[0] != -1 && genotypes[1] != -1 ) {
				++table( genotypes[0][i], genotypes[1][i] ) ;
			}
		}
		
		
	}

	void end_processing_snps() {}
	
	typedef boost::function< void( genfile::SNPIdentifyingData const& source, genfile::SNPIdentifyingData const& target, double, double, double ) > ResultCallback ;
	
private:
	
	genfile::SNPDataSource::UniquePtr m_source ;
	double const m_threshhold ;
	boost::function< void( genfile::SNPIdentifyingData const& source, genfile::SNPIdentifyingData const& target, double, double, double ) > m_result_callback ;
} ;

#endif
