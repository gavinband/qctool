#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/function.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/vcf/get_set.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "SNPSummaryComponent.hpp"

void SNPSummaryComputationManager::add_computation( std::string const& name, SNPSummaryComputation::UniquePtr computation ) {
	m_computations.insert( name, computation ) ;
}

void SNPSummaryComputationManager::add_result_callback( ResultCallback callback ) {
	m_result_signal.connect( callback ) ;
}

void SNPSummaryComputationManager::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	m_snp_index = 0 ;
}

void SNPSummaryComputationManager::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	{
		genfile::vcf::MatrixSetter< SNPSummaryComputation::Genotypes > setter( m_genotypes ) ;
		data_reader.get( "genotypes", setter ) ;
	}
	Computations::const_iterator i = m_computations.begin(), end_i = m_computations.end() ;
	for( ; i != end_i; ++i ) {
		i->second->operator()(
			snp,
			m_genotypes,
			boost::bind(
				boost::ref( m_result_signal ),
				m_snp_index,
				snp,
				i->first,
				_1,
				_2
			)
		) ;
	}
	++m_snp_index ;
}

void SNPSummaryComputationManager::end_processing_snps() {}

void SNPSummaryComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "SNP computation options" ) ;
	options[ "-snp-stats" ]
		.set_description( "Calculate and output per-SNP statistics.  This implies that no SNP filtering options are used." ) ;
    options[ "-snp-stats-file" ]
        .set_description( 	"Override the auto-generated path(s) of the snp-stats file to use when outputting snp-wise statistics.  "
							"(By default, the paths are formed by adding \".snp-stats\" to the input gen filename(s).)  "
							"The '#' character can also be used here to specify one output file per chromosome." )
        .set_takes_values_per_use(1)
		.set_maximum_multiplicity(1)
		.set_default_value( "qctool.snp-stats" ) ;

	options[ "-snp-stats-columns" ]
        .set_description( "Comma-seperated list of extra columns to output in the snp-wise statistics file.  "
 						"The standard columns are: "
						"SNPID, RSID, position, minor_allele, major_allele, MAF, HWE, missing, information."
						" Your choices here are old_information, jonathans_information, mach_r2, and entropy." )
		.set_takes_single_value()
		.set_default_value( "alleles,HWE,missingness,information" ) ;
}

SNPSummaryComponent::SNPSummaryComponent( appcontext::OptionProcessor const& options ):
	m_options( options )
{}

namespace {
	struct FileOutputter: public boost::noncopyable {
		typedef std::auto_ptr< FileOutputter > UniquePtr ;
		typedef boost::shared_ptr< FileOutputter > SharedPtr ;
		
		static UniquePtr create( std::string const& filename ) { return UniquePtr( new FileOutputter( filename ) ) ; }

		FileOutputter( std::string const& filename ):
			m_filename( filename ),
			m_sink( statfile::BuiltInTypeStatSink::open( filename ))
		{
			(*m_sink) | "SNPID" | "rsid" | "chromosome" | "position" | "alleleA" | "alleleB" | "computation_name" | "value_name" | "value" ;
		}

		void operator()(
			std::size_t index,
			genfile::SNPIdentifyingData const& snp,
			std::string const& computation_name,
			std::string const& value_name,
			genfile::VariantEntry const& value
		) {
			(*m_sink)
				<< snp.get_SNPID()
				<< snp.get_rsid()
				<< std::string( snp.get_position().chromosome() )
				<< snp.get_position().position()
				<< snp.get_first_allele()
				<< snp.get_second_allele()
				<< computation_name
				<< value_name
				<< value
				<< statfile::end_row() ;
			;
		}

	private:
		std::string const m_filename ;
		statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
	} ;
}

genfile::SNPDataSourceProcessor::Callback::UniquePtr SNPSummaryComponent::create() const {
	genfile::SNPDataSourceProcessor::Callback::UniquePtr result( create_manager().release() ) ;
	return result ;
}
	

SNPSummaryComputationManager::UniquePtr SNPSummaryComponent::create_manager() const {
	SNPSummaryComputationManager::UniquePtr manager( new SNPSummaryComputationManager() ) ;
	std::vector< std::string > elts = genfile::string_utils::split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "-snp-stats-columns" ), ",", " \t" ) ;
	foreach( std::string const& elt, elts ) {
		manager->add_computation( elt, SNPSummaryComputation::create( elt )) ;
	}
	manager->add_result_callback(
		boost::bind(
			&FileOutputter::operator(),
			FileOutputter::SharedPtr( FileOutputter::create( m_options.get_value< std::string >( "-snp-stats-file" )) ),
			_1, _2, _3, _4, _5
		)
	) ;
	return manager ;
}

SNPSummaryComputation::UniquePtr SNPSummaryComponent::create_computation( std::string const& name ) const {
	return SNPSummaryComputation::UniquePtr( SNPSummaryComputation::create( name )) ;
}
