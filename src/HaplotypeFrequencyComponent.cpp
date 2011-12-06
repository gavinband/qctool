#include <string>
#include <boost/function.hpp>
#include <Eigen/Core>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/RandomAccessSNPDataSource.hpp"
#include "genfile/SNPDataSourceChain.hpp"
#include "genfile/wildcard.hpp"
#include "genfile/vcf/get_set.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "HaplotypeFrequencyComponent.hpp"
#include "integration/NewtonRaphson.hpp"
#include "integration/Derivative.hpp"

namespace impl {
	struct HaplotypeFrequencyFileOutputter ;
}

struct HaplotypeFrequencyLogLikelihood ;

void HaplotypeFrequencyComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "LD computation options" ) ;
	options[ "-compute-ld-with" ]
		.set_description( "Compute LD pairwise metrics between the main dataset and SNPs." )
		.set_takes_single_value() ;
	options[ "-compute-ld-file" ]
		.set_description( "File in which to place computation of pairwise SNP LD measures." )
		.set_takes_single_value()
		.set_default_value( "ld.txt" ) ;
	options.option_implies_option( "-compute-ld-file", "-compute-ld-with" ) ;
}

namespace impl {
	struct HaplotypeFrequencyFileOutputter {
	public:
		typedef std::auto_ptr< HaplotypeFrequencyFileOutputter > UniquePtr ;
		typedef boost::shared_ptr< HaplotypeFrequencyFileOutputter > SharedPtr ;
		static UniquePtr create( std::string const& filename ) {
			return UniquePtr( new HaplotypeFrequencyFileOutputter( filename )) ;
		}
		static SharedPtr create_shared( std::string const& filename ) {
			return SharedPtr( new HaplotypeFrequencyFileOutputter( filename )) ;
		}
		HaplotypeFrequencyFileOutputter( std::string filename ): m_filename( filename ) {}
		
		void operator()(
			genfile::SNPIdentifyingData const& source_snp,
			genfile::SNPIdentifyingData const& target_snp,
			Eigen::VectorXd const& measures
		) {
			if( !m_sink.get() ) {
				m_sink = statfile::BuiltInTypeStatSink::open( m_filename ) ;
				(*m_sink) | "SNPID1" | "rsid1" | "chromosome1" | "position1" | "SNPID2" | "rsid2" | "chromsome2" | "position2" | "pi_00" | "pi_01" | "pi_10" | "pi_11" ;
			}
			(*m_sink)
				<< source_snp.get_SNPID() << source_snp.get_rsid() << source_snp.get_position().chromosome() << source_snp.get_position().position()
				<< target_snp.get_SNPID() << target_snp.get_rsid() << target_snp.get_position().chromosome() << target_snp.get_position().position()
				<< ( 1.0 - measures.sum() ) << measures(0) << measures(1) << measures(2)
				<< statfile::end_row() ;
		}
		
	private:
		std::string const m_filename ;
		statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
	} ;
}

HaplotypeFrequencyComponent::UniquePtr HaplotypeFrequencyComponent::create( appcontext::OptionProcessor const& options ) {
	HaplotypeFrequencyComponent::UniquePtr result ;
	result.reset(
		new HaplotypeFrequencyComponent(
			genfile::SNPDataSource::create_chain(
				genfile::wildcard::find_files_by_chromosome(
					options.get< std::string >( "-compute-ld-with" ),
					genfile::wildcard::eALL_CHROMOSOMES
				)
			)
		)
	) ;
	
	impl::HaplotypeFrequencyFileOutputter::SharedPtr outputter = impl::HaplotypeFrequencyFileOutputter::create_shared(
		options.get_value< std::string >( "-compute-ld-file" )
	) ;
	
	result->send_results_to( boost::bind( &impl::HaplotypeFrequencyFileOutputter::operator(), outputter, _1, _2, _3 ) ) ;
	return result ;
}

HaplotypeFrequencyComponent::HaplotypeFrequencyComponent(
	genfile::SNPDataSource::UniquePtr source
):
	m_source( source ),
	m_threshhold( 0.9 )
{}

void HaplotypeFrequencyComponent::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	assert( m_source->number_of_samples() == number_of_samples ) ;
}

void HaplotypeFrequencyComponent::processed_snp( genfile::SNPIdentifyingData const& target_snp, genfile::VariantDataReader& target_data_reader ) {
	genfile::SNPIdentifyingData source_snp ;
	genfile::SingleSNPGenotypeProbabilities source_probs, target_probs ;
	target_data_reader.get( "genotypes", target_probs ) ;
	m_source->reset_to_start() ;
	while( m_source->get_snp_identifying_data( source_snp )) {
		genfile::VariantDataReader::UniquePtr source_data_reader = m_source->read_variant_data() ;
		compute_ld_measures( source_snp, *source_data_reader, target_snp, target_data_reader ) ;
	}
}

void HaplotypeFrequencyComponent::compute_ld_measures(
	genfile::SNPIdentifyingData const& source_snp,
	genfile::VariantDataReader& source_data_reader,
	genfile::SNPIdentifyingData const& target_snp,
	genfile::VariantDataReader& target_data_reader
) {
	std::vector< std::vector< int > > genotypes( 2 ) ;
	genfile::vcf::GenotypeSetter< std::vector< int > > source_getter( genotypes[0], m_threshhold ) ;
	genfile::vcf::GenotypeSetter< std::vector< int > > target_getter( genotypes[1], m_threshhold ) ;
	source_data_reader.get( "genotypes", source_getter ) ;
	target_data_reader.get( "genotypes", target_getter ) ;
	assert( genotypes[0].size() == m_source->number_of_samples() ) ;
	assert( genotypes[0].size() == genotypes[1].size() ) ;
	compute_ld_measures(
		source_snp,
		target_snp,
		genotypes
	) ;
}

void HaplotypeFrequencyComponent::compute_ld_measures(
	genfile::SNPIdentifyingData const& source_snp,
	genfile::SNPIdentifyingData const& target_snp,
	std::vector< std::vector< int > > const& genotypes
) {
	// Construct table of genotypes at each SNP.
	Eigen::Matrix3d table = Eigen::Matrix3d::Zero() ;
	for( std::size_t i = 0; i < m_source->number_of_samples(); ++i ) {
		if( genotypes[0][i] != -1 && genotypes[1][i] != -1 ) {
			++table( genotypes[0][i], genotypes[1][i] ) ;
		}
	}
	
	HaplotypeFrequencyLogLikelihood ll( table ) ;
	typedef HaplotypeFrequencyLogLikelihood::Vector Vector ;
	Vector params( 3 ) ;
	params << 0.2, 0.5, 0.2 ;
	
	using integration::derivative ;
	integration::Derivative< HaplotypeFrequencyLogLikelihood > Dll = derivative( ll ) ;
	
	params = integration::find_root_by_newton_raphson(
		Dll,
		params,
		0.00001
	) ;

	m_result_signal( source_snp, target_snp, params ) ;
}

void HaplotypeFrequencyComponent::end_processing_snps() {}

void HaplotypeFrequencyComponent::send_results_to( ResultCallback callback ) {
	m_result_signal.connect( callback ) ;
}
