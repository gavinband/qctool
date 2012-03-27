#include <string>
#include <memory>
#include <set>
#include <map>
#include <utility>
#include <boost/bind.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/function.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/regex.hpp>
#include "appcontext/CmdLineOptionProcessor.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "appcontext/ProgramFlow.hpp"
#include "genfile/FileUtils.hpp"
#include "genfile/utility.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include <Eigen/Core>

namespace globals {
	std::string const program_name = "amet" ;
	std::string const program_version = "0.1" ;
}

struct FrequentistGenomeWideAssociationResults: public boost::noncopyable {
	typedef std::auto_ptr< FrequentistGenomeWideAssociationResults > UniquePtr ;
	virtual ~FrequentistGenomeWideAssociationResults() {}
	typedef boost::function< void ( std::size_t i, genfile::SNPIdentifyingData const& snp ) > SNPCallback ;
	virtual void get_SNPs( SNPCallback ) const = 0 ;
	virtual void get_results(
		std::size_t snp_i,
		boost::function< void ( double, double, double, double ) > get_counts ,
		boost::function< void ( genfile::VariantEntry ) > get_allele_frequency,
		boost::function< void ( genfile::VariantEntry ) > get_beta,
		boost::function< void ( genfile::VariantEntry ) > get_se
	) = 0 ;
	virtual std::string get_summary( std::string const& prefix = "", std::size_t target_column = 80 ) const = 0;
} ;


struct SNPTESTResults: public FrequentistGenomeWideAssociationResults {
	typedef boost::function< void ( std::size_t, boost::optional< std::size_t > ) > ProgressCallback ;
	SNPTESTResults( std::vector< genfile::wildcard::FilenameMatch > const& filenames, ProgressCallback progress_callback = ProgressCallback() )
	{
		setup( filenames, progress_callback ) ;
	}

	void get_SNPs( SNPCallback callback ) const {
		for( std::size_t i = 0; i < m_snps.size(); ++i ) {
			callback( i, m_snps[i] ) ;
		}
	}

	virtual void get_results(
		std::size_t snp_i,
		boost::function< void ( double, double, double, double ) > get_counts ,
		boost::function< void ( genfile::VariantEntry ) > get_allele_frequency,
		boost::function< void ( genfile::VariantEntry ) > get_beta,
		boost::function< void ( genfile::VariantEntry ) > get_se
	) {
		assert( snp_i < m_snps.size() ) ;
		get_counts( m_sample_counts( snp_i, 0 ), m_sample_counts( snp_i, 1 ), m_sample_counts( snp_i, 2 ), m_sample_counts( snp_i, 3 ) ) ;
		get_allele_frequency( ( m_sample_counts( snp_i, 1 ) + 2 * m_sample_counts( snp_i, 2 ) ) / ( 2.0 * m_sample_counts.row( snp_i ).head( 3 ).sum() )) ;
		get_beta( m_betas( snp_i ) ) ;
		get_se( m_ses( snp_i ) ) ;
	}

	std::string get_summary( std::string const& prefix, std::size_t target_column ) const {
		using genfile::string_utils::to_string ;
		std::string result = prefix + "SNPTESTResults object holding results for " + to_string( m_snps.size() ) + " SNPs.  The first few are:\n" ;
		for( std::size_t i = 0; i < std::min( std::size_t( 5 ), m_snps.size() ); ++i ) {
			result += to_string( m_snps[i] ) + ": beta=" + to_string( m_betas[i] ) + ", se=" + to_string( m_ses[i] ) + ", pvalue=" + to_string( m_pvalues[i] ) + "\n" ;
		}
		return result ;
	}

	private:
		std::vector< genfile::SNPIdentifyingData > m_snps ;
		Eigen::VectorXd m_betas ;
		Eigen::VectorXd m_ses ;
		Eigen::VectorXd m_pvalues ;
		Eigen::VectorXd m_info ;
		Eigen::MatrixXd m_sample_counts ;
		
		void setup( std::vector< genfile::wildcard::FilenameMatch > const& filenames, ProgressCallback progress_callback = ProgressCallback() ) {
			for( std::size_t i = 0; i < filenames.size(); ++i ) {
				setup( filenames[i] ) ;
				progress_callback( i+1, filenames.size() ) ;
			}
		}
		
		void setup( genfile::wildcard::FilenameMatch const& filename ) {
			statfile::BuiltInTypeStatSource::UniquePtr source( statfile::BuiltInTypeStatSource::open( filename.filename() )) ;
			std::map< std::string, std::size_t > columns_by_name ;
			std::map< std::size_t, std::string > columns_by_index ;
			columns_by_name[ "_beta_1" ] = 0 ;
			columns_by_name[ "_se_1" ] = 0 ;
			columns_by_name[ "_pvalue" ] = 0 ;
			columns_by_name[ "info" ] = 0 ;
			columns_by_name[ "all_AA" ] = 0 ;
			columns_by_name[ "all_AB" ] = 0 ;
			columns_by_name[ "all_BB" ] = 0 ;
			columns_by_name[ "all_NULL" ] = 0 ;
			for( std::size_t i = 0; i < source->number_of_columns(); ++i ) {
				std::string name = source->name_of_column( i ) ;
				for( std::map< std::string, std::size_t >::iterator j = columns_by_name.begin(); j != columns_by_name.end(); ++j ) {
					if( name.size() >= j->first.size() && name.compare( name.size() - j->first.size(), j->first.size(), j->first ) == 0 ) {
						j->second = i ;
						columns_by_index[i] = j->first ;
						
						std::cerr << "Column \"" << j->first << "\" found at position " << i << ".\n" ;
					}
				}
			}
		
			m_snps.resize( source->number_of_rows() ) ;
			m_betas.resize( source->number_of_rows() ) ;
			m_pvalues.resize( source->number_of_rows() ) ;
			m_ses.resize( source->number_of_rows() ) ;
			m_info.resize( source->number_of_rows() ) ;
			m_sample_counts.resize( source->number_of_rows(), 4 ) ;

			using genfile::string_utils::to_repr ;
			
			genfile::SNPIdentifyingData snp ;
			for(
				std::size_t index = 0 ;
				(*source) >> snp.SNPID() >> snp.rsid() >> snp.position().chromosome() >> snp.position().position() >> snp.first_allele() >> snp.second_allele();
				++index, (*source) >> statfile::ignore_all()
			) {
				m_snps[ index ] = snp ;
		
				std::map< std::size_t, std::string >::const_iterator
					i = columns_by_index.begin(),
					end_i = columns_by_index.end() ;
				for( ; i != end_i; ++i ) {
					std::string value ;
					(*source) >> statfile::ignore( i->first - source->current_column() ) >> value ;
					
					if( i->second == "_beta_1" ) {
						m_betas( index ) = to_repr< double >( value ) ;
					}
					else if( i->second == "_se_1" ) {
						m_ses( index ) = to_repr< double >( value ) ;
					}
					else if( i->second == "_pvalue" ) {
						m_pvalues( index ) = to_repr< double >( value ) ;
					}
					else if( i->second == "info" ) {
						m_info( index ) = to_repr< double >( value ) ;
					}
					else if( i->second == "all_AA" ) {
						m_sample_counts( index, 0 ) = to_repr< double >( value ) ;
					}
					else if( i->second == "all_AB" ) {
						m_sample_counts( index, 1 ) = to_repr< double >( value ) ;
					}
					else if( i->second == "all_BB" ) {
						m_sample_counts( index, 2 ) = to_repr< double >( value ) ;
					}
					else if( i->second == "all_NULL" ) {
						m_sample_counts( index, 3 ) = to_repr< double >( value ) ;
					}
				}
			}
		}
} ;

struct AmetOptions: public appcontext::CmdLineOptionProcessor {
	std::string get_program_name() const { return globals::program_name ; }

	void declare_options( appcontext::OptionProcessor& options ) {
		options.set_help_option( "-help" ) ;

		options.declare_group( "File handling options" ) ;
		options[ "-snptest" ]
			.set_description( "Specify the path of a file containing SNPTEST results to load." )
			.set_is_required()
			.set_takes_values( 1 )
			.set_minimum_multiplicity( 1 )
			.set_maximum_multiplicity( 100 )
		;
		
		options[ "-snp-match-fields" ]
			.set_description( "Set the fields used to describe SNPs as equal." )
			.set_takes_single_value()
			.set_default_value( "position,alleles" ) ;

		options[ "-log" ]
			.set_description( "Specify the path of the log file." )
			.set_takes_single_value()
			.set_default_value( "overrep.log" ) ;
	}
} ;


struct AmetApplication: public appcontext::ApplicationContext {
public:
	AmetApplication( int argc, char **argv ):
		appcontext::ApplicationContext(
			globals::program_name,
			globals::program_version,
			std::auto_ptr< appcontext::OptionProcessor >( new AmetOptions ),
			argc,
			argv,
			"-log"
		),
		m_snps( genfile::SNPIdentifyingData::CompareFields( options().get< std::string > ( "-snp-match-fields" )) )
	{}
	
	void process() {
		unsafe_process() ;
	}
	
private:
	boost::ptr_vector< FrequentistGenomeWideAssociationResults > m_cohorts ;
	typedef std::map< genfile::SNPIdentifyingData, std::vector< boost::optional< std::size_t > >, genfile::SNPIdentifyingData::CompareFields > SnpMap ;
	SnpMap m_snps ;
	typedef std::map< std::vector< bool >, std::size_t > CategoryCounts ;
	CategoryCounts m_category_counts ;
private:
	
	void unsafe_process() {
		load_data() ;
		link_data() ;
		categorise_snps() ;
		summarise() ;
	}
	
	void summarise() {
		get_ui_context().logger() << "================================================\n" ;
		get_ui_context().logger() << "Cohort summary:\n" ;
		for( std::size_t i = 0; i < m_cohorts.size(); ++i ) {
			get_ui_context().logger() << m_cohorts[i].get_summary( " - " ) ;
		}
		
		get_ui_context().logger() << "\n================================================\n" ;
		get_ui_context().logger() << "SNP Categories:\n" ;
		for( CategoryCounts::const_iterator i = m_category_counts.begin(); i != m_category_counts.end(); ++i ) {
			for( std::size_t j = 0; j < m_cohorts.size(); ++j ) {	
				get_ui_context().logger() << " " << i->first[j] ;
			}
			get_ui_context().logger() << ": " << i->second << "\n" ;
		}
		get_ui_context().logger() << " TOTAL: " << m_snps.size() << "\n" ;
		get_ui_context().logger() << "================================================\n" ;
	}

	void add_SNP_callback( std::size_t cohort_i, std::size_t snp_i, genfile::SNPIdentifyingData const& snp ) {
		std::vector< boost::optional< std::size_t > >& snp_indices = m_snps[ snp ] ;
		snp_indices.resize( m_cohorts.size() ) ;
		snp_indices[ cohort_i ] = snp_i ;
	}

	void link_data() {
		UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Linking SNPs" ) ;
		progress_context( 0, m_cohorts.size() ) ;
		for( std::size_t cohort_i = 0; cohort_i < m_cohorts.size(); ++cohort_i ) {
			m_cohorts[ cohort_i].get_SNPs(
				boost::bind(
					&AmetApplication::add_SNP_callback,
					this,
					cohort_i,
					_1,
					_2
				)
			) ;
			progress_context( cohort_i+1, m_cohorts.size() ) ;
		}
	}
	
	void categorise_snps() {
		SnpMap::const_iterator snp_i = m_snps.begin(), end_i = m_snps.end() ;
		for( ; snp_i != end_i; ++snp_i ) {
			std::vector< boost::optional< std::size_t > > const& indices = snp_i->second ;
			std::vector< bool > indicator( indices.size(), false ) ;
			for( std::size_t i = 0; i < indices.size(); ++i ) {
				indicator[i] = bool( indices[i] ) ;
			}
			++m_category_counts[ indicator ] ;
		}
	}

	void load_data() {
		std::vector< std::string > cohort_files = options().get_values< std::string >( "-snptest" ) ;
		for( std::size_t i = 0; i < cohort_files.size(); ++i ) {
			UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading SNPTEST results \"" + cohort_files[i] + "\"" ) ;
			FrequentistGenomeWideAssociationResults::UniquePtr results( new SNPTESTResults( genfile::wildcard::find_files_by_chromosome( cohort_files[i] ), progress_context ) ) ;
			m_cohorts.push_back( results ) ;
		}
	}
} ;

int main( int argc, char **argv ) {
	try {
		AmetApplication app( argc, argv ) ;	
		app.process() ;
	}
	catch( appcontext::HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}
