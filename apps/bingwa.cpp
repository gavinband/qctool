
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <sstream>
#include <iomanip>
#include <memory>
#include <set>
#include <map>
#include <utility>
#include <boost/bimap.hpp>
#include <boost/bind.hpp>
#include <boost/signals2/signal.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/function.hpp>
#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/regex.hpp>
#include <boost/math/distributions/normal.hpp>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include "appcontext/CmdLineOptionProcessor.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "appcontext/ProgramFlow.hpp"
#include "genfile/FileUtils.hpp"
#include "genfile/utility.hpp"
#include "genfile/string_utils/string_utils.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPIdentifyingData2.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/CommonSNPFilter.hpp"
#include "genfile/VariantEntry.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include "statfile/SNPDataSourceAdapter.hpp"
#include "qcdb/Storage.hpp"
#include "qcdb/FlatFileOutputter.hpp"
#include "qcdb/FlatTableDBOutputter.hpp"
#include "FrequentistGenomeWideAssociationResults.hpp"
#include "EffectParameterNamePack.hpp"
#include "SNPTESTResults.hpp"
#include "MMMResults.hpp"

// #define DEBUG_BINGWA 1

namespace globals {
	std::string const program_name = "bingwa" ;
	std::string const program_version = "0.2" ;
}

namespace {
	double const NA = std::numeric_limits< double >::quiet_NaN() ;
	double const Infinity = std::numeric_limits< double >::infinity() ;
}


struct BingwaOptions: public appcontext::CmdLineOptionProcessor {
	std::string get_program_name() const { return globals::program_name ; }

	void declare_options( appcontext::OptionProcessor& options ) {
		options.set_help_option( "-help" ) ;

		{
			options.declare_group( "Input data options" ) ;
			options[ "-data" ]
				.set_description( "Specify the path of a file containing SNPTEST or MMM results to load." )
				.set_is_required()
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 1 )
				.set_maximum_multiplicity( 100 )
			;
			options[ "-extra-columns" ]
				.set_description( "Specify extra columns in input files whose values will be considered as variables to be reported in the output."
				 	" Currently these must be columns of numerical data.  (A single wildcard character * at the start or end of the column name may"
					" be used to match any initial or terminal sequence of characters; but take care to escape this from the shell.)" )
					.set_takes_values_until_next_option()
					.set_minimum_multiplicity( 0 )
					.set_maximum_multiplicity( 100 )
				;
			
			options[ "-snp-match-fields" ]
				.set_description( "Use this option to specify a comma-separated list of SNP-identifying fields that should be used "
					"to match SNPs between cohorts.  Possible fields "
					"are \"position\", \"alleles\", \"rsid\", or \"snpid\"; you must always specify \"position\" as the first entry." )
				.set_takes_single_value()
				.set_default_value( "position,alleles" ) ;
		}
		{
			options.declare_group( "SNP inclusion / exclusion options" ) ;
			options[ "-excl-snpids" ]
				.set_description( "Exclude all SNPs whose SNPID is in the given file(s) from the analysis.")
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-excl-rsids" ]
				.set_description( "Exclude all SNPs whose RSID is in the given file(s) from the analysis.")
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-incl-snpids" ]
				.set_description( "Exclude all SNPs whose SNPID is not in the given file(s) from the analysis.")
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-incl-rsids" ]
				.set_description( "Exclude all SNPs whose RSID is not in the given file(s) from the analysis.")
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-excl-positions" ]
				.set_description( "Exclude all SNPs whose position is in the given file(s) from the analysis. "
					"Positions should be in the form [chromosome]:[position] and separated by whitespace." )
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-incl-positions" ]
				.set_description( "Exclude all SNPs whose position is not in the given file(s) from the analysis. "
					"Positions should be in the form [chromosome]:[position] and separated by whitespace." )
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-excl-snps-matching" ]
				.set_description( "Filter out snps whose rsid or SNPID matches the given value. "
					"The value should be a string which can contain a % wildcard character (which matches any substring). "
					"Optionally, prefix the argument with snpid~ or rsid~ to only match against the SNPID or rsid fields." )
				.set_takes_single_value() ;
			options[ "-incl-snps-matching" ]
				.set_description( "Filter out snps whose rsid or SNPID does not match the given value. "
					"The value should be a string which can contain a % wildcard character (which matches any substring). "
					"Optionally, prefix the argument with snpid~ or rsid~ to only match against the SNPID or rsid fields." )
				.set_takes_single_value() ;
			options[ "-incl-range" ]
				.set_description( "Specify a range of SNPs (or comma-separated list of ranges of SNPs) to operate on. "
					"Each range should be in the format CC:xxxx-yyyy where CC is the chromosome and xxxx and yyyy are the "
					"start and end coordinates, or just xxxx-yyyy which matches that range from all chromosomes. "
					"You can also omit either of xxxx or yyyy to get all SNPs from the start or to the end of a chromosome." )
				.set_takes_single_value() ;
			options[ "-excl-range" ]
				.set_description( "Specify a range of SNPs (or comma-separated list of ranges of SNPs) to exclude from operation. "
					"Each range should be in the format CC:xxxx-yyyy where CC is the chromosome and xxxx and yyyy are the "
					"start and end coordinates, or just xxxx-yyyy which matches that range from all chromosomes. "
					"You can also omit either of xxxx or yyyy to get all SNPs from the start or to the end of a chromosome." )
				.set_takes_single_value() ;
				
			options[ "-excl-snpids-per-cohort" ]
				.set_description( "Exclude all SNPs whose SNPID is in the given file(s) from the corresponding cohort. "
					"The i-th value should be a comma-separated list of files containing rsids to exclude from the i-th cohort." )
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-excl-rsids-per-cohort" ]
				.set_description( "Exclude all SNPs whose RSID is in the given file(s) from the corresponding cohort. "
					"The i-th value should be a comma-separated list of files containing rsids to exclude from the i-th cohort." )
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-incl-snpids-per-cohort" ]
				.set_description( "Exclude all SNPs whose SNPID is not in the given file(s) from the corresponding cohort."
					"The i-th value should be a comma-separated list of files containing rsids to exclude from the i-th cohort." )
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-incl-rsids-per-cohort" ]
				.set_description( "Exclude all SNPs whose RSID is not in the given file(s) from the corresponding cohort."
					"The i-th value should be a comma-separated list of files containing rsids to exclude from the i-th cohort." )
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-excl-snps-per-cohort" ]
				.set_description( "Exclude all SNPs in the given file(s) from the corresponding cohort. "
					"The i-th value should be a comma-separated list of files containing snps to exclude from the i-th cohort, in the format output by gen-grep." )
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-incl-snps-per-cohort" ]
				.set_description( "Exclude all SNPs not in the given file(s) from the corresponding cohort."
					"The i-th value should be a comma-separated list of files containing rsids to exclude from the i-th cohort, in the format output by gen-grep." )
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
		
			options[ "-min-info" ]
				.set_description( "Treat SNPs with info less than the given threshhold as missing." )
				.set_takes_values( 1 ) ;

			options[ "-min-maf" ]
				.set_description( "Treat SNPs with maf (in controls) less than the given threshhold as missing." )
				.set_takes_values( 1 ) ;
		}

		{
			options.declare_group( "Options for adjusting variants" ) ;
			options[ "-flip-alleles" ]
				.set_description( "Specify that, where alleles do not match between cohorts, bingwa will try to match them by flipping." )
			;
			options[ "-assume-chromosome" ]
				.set_description( "Specify that bingwa should treat all SNPs with missing chromosome as having the given one." )
				.set_takes_single_value() ;
		}

		{
			options.declare_group( "Options affecting output" ) ;

			options[ "-o" ]
				.set_description( "Specify the path to the output file." )
				.set_is_required()
				.set_takes_single_value() ;

			options[ "-flat-file" ]
				.set_description( "Specify the output file should be a flat file, not a db." ) ;
			options[ "-noindex" ]
				.set_description( "Specify that " + globals::program_name + " should not create large indices on database tables when finalising storage."
					" Indices are usually desired.  However, when running very large jobs parallelised across subsets, it is generally faster to"
					" use -noindex and create indices manually when all jobs have completed." ) ;
			options[ "-analysis-name" ]
				.set_description( "Specify a name for the current analysis." )
				.set_takes_single_value()
				.set_default_value( "bingwa analysis" ) ;
			options[ "-analysis-chunk" ]
				.set_description( "Specify a name denoting the current genomic region or chunk on which this is run.  This is intended for use in parallel environments." )
				.set_takes_single_value()
				.set_default_value( genfile::MissingValue() ) ;
			options[ "-cohort-names" ]
				.set_description( "Specify a name to label results from this analysis with" )
				.set_takes_values_until_next_option() ;
			options[ "-table-name" ]
				.set_description( "Specify a name for the table to use." )
				.set_takes_single_value() ;

			options[ "-log" ]
				.set_description( "Specify the path of a log file; all screen output will be copied to the file." )
				.set_takes_single_value() ;
		}

		{
			options.declare_group( "Analysis options" ) ;
		
			options[ "-no-meta-analysis" ]
				.set_description( "Don't do a fixed effect meta-analysis.  Instead, just match up SNPs and store per-cohort values." ) ;

			options[ "-define-sd-set" ]
				.set_description( "Define set(s) of standard deviations to use in prior specifications. "
					"The value should be a whitespace-separated list of elements of the form "
					"\"[name of set]=[sd1],[sd2],[sd3]...\"." )
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 0 )
				.set_maximum_multiplicity( 1 ) ;

			options[ "-define-rho-set" ]
				.set_description( "Define set(s) of correlations to use in prior specifications. "
					"The value should be a whitespace-separated list of elements of the form "
					"\"[name of set]=[value],[value],[value]...\"." )
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 0 )
				.set_maximum_multiplicity( 1 ) ;
			
			options[ "-simple-prior" ]
				.set_description( "Specify a model of the form \"rho=[r]/sd=[s]\" giving a correlation matrix of the form\n"
					"      [ 1 r r  .. ]\n"
					"      [ r 1 r  .. ]\n"
					"s^2 * [ r r 1  .. ]\n"
					"      [       .   ]\n"
					"      [         . ]\n"
				)
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 0 )
				.set_maximum_multiplicity( 1 ) ;

			options[ "-simple-prior-name" ]
				.set_description( "Specify a name for each simple prior." )
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 0 )
				.set_maximum_multiplicity( 1 ) ;
			
			options[ "-complex-prior" ]
				.set_description( "Specify the upper triangle (including diagonal) of the prior correlation matrix for bayesian analysis. "
					"For example, the value \"1,0.5,1\" specifies the matrix\n"
					"   [ 1    0.5 ]\n"
					"   [ 0.5  1   ]."
				)
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 0 )
				.set_maximum_multiplicity( 1 ) ;
			
			options[ "-complex-prior-name" ]
				.set_description( "Specify the name of the complex models.  This option takes the same number of values as "
					"are given to -complex-prior option."
				)
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 0 )
				.set_maximum_multiplicity( 1 )
			;

			options[ "-define-group" ]
				.set_description( "Specify a group of cohorts for use with -group-prior and -group-specific prior."
					" The format for each group is [group name]=[cohort name 1],[cohort name 2],..." )
					.set_takes_values_until_next_option()
					.set_minimum_multiplicity( 0 )
					.set_maximum_multiplicity( 1 )
				;
			
			options[ "-group-prior" ]
				.set_description( "Specify a prior made up of blocks per groups, in the format"
					" [group1]:rho=[r1]/sd=[s1],[group2]:rho=[r2]/sd=[s2]...:rho=[r]" )
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 0 )
				.set_maximum_multiplicity( 1 )
			;

			options[ "-group-prior-name" ]
				.set_description( "Specify the name of the group priors." )
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 0 )
				.set_maximum_multiplicity( 1 )
			;

			options[ "-group-specific-prior" ]
				.set_description( "Specify a group-specific prior in the format [cohorts]:sd=[value]/rho=[value]" )
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 0 )
				.set_maximum_multiplicity( 1 )
			;
			
			options[ "-group-specific-prior-name" ]
				.set_description( "Specify the name of the group-specific priors." )
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 0 )
				.set_maximum_multiplicity( 1 )
			;
			
			options.option_excludes_option( "-no-meta-analysis", "-simple-prior" ) ;
			options.option_excludes_option( "-no-meta-analysis", "-complex-prior" ) ;
			options.option_excludes_option( "-no-meta-analysis", "-group-specific-prior" ) ;

			options.option_implies_option( "-simple-prior-name", "-simple-prior" ) ;

			options.option_implies_option( "-complex-prior", "-complex-prior-name" ) ;
			options.option_implies_option( "-complex-prior-name", "-complex-prior" ) ;
//			options.option_implies_option( "-group-prior", "-define-group" ) ;
//			options.option_implies_option( "-group-specific-prior", "-define-group" ) ;
			options.option_implies_option( "-group-specific-prior-name", "-group-specific-prior" ) ;
		}
	}
} ;

FrequentistGenomeWideAssociationResults::UniquePtr FrequentistGenomeWideAssociationResults::create(
	std::vector< genfile::wildcard::FilenameMatch > const& filenames,
	boost::optional< std::string > const& effect_size_column_regex,
	std::vector< std::string > const& columns,
	genfile::SNPIdentifyingDataTest::UniquePtr test,
	boost::optional< genfile::Chromosome > chromosome_hint,
	SNPResultCallback result_callback,
	ProgressCallback progress_callback
) {
	// peek at the file to determine file type
	
	std::string type = "unknown" ;

	if( filenames.size() > 0 ) {
		std::auto_ptr< std::istream > file = genfile::open_text_file_for_input( filenames[0].filename() ) ;
		std::string line ;
		std::getline( *file, line ) ;
		std::vector< std::string > const elts = genfile::string_utils::split( line, " \t," ) ;
		if( elts.size() > 4
			&& elts[0] == "chr"
			&& elts[1] == "snp_id1"
			&& elts[2] == "snp_id2"
			&& elts[3] == "pos"
			&& elts[4] == "allele_0"
		) {
			type = "mmm" ; // Matti's mixed model, http://www.well.ox.ac.uk/~mpirinen/
		}
		else if(
			elts.size() > 5
			&& elts[0] == "snp_id1"
			&& elts[1] == "snp_id2"
			&& elts[2] == "pos"
			&& elts[3] == "allele_0"
			&& elts[4] == "allele_1"
			&& elts[5] == "n_included"
		) {
			type = "mmm" ; // Matti's mixed model, http://www.well.ox.ac.uk/~mpirinen/
		}
		else if( line.substr( 0, 40 ) == "id rsid chromosome pos allele_A allele_B" ) {
				type = "snptest" ;
		}
	}

	FlatFileFrequentistGenomeWideAssociationResults::UniquePtr result ;
	if( type == "snptest" || type == "unknown" ) {
		result.reset( new SNPTESTResults( test, chromosome_hint ) ) ; 
	}
	else if( type == "mmm" ) {
		result.reset( new MMMResults( test, chromosome_hint ) ) ;
	}
	
	if( effect_size_column_regex ) {
		result->set_effect_size_column_regex( effect_size_column_regex.get() ) ;
	}

	BOOST_FOREACH( std::string const& column, columns ) {
		result->add_variable( column ) ;
	}

	result->add_data(
		filenames,
		result_callback,
		progress_callback
	) ;

	return FrequentistGenomeWideAssociationResults::UniquePtr( result.release() ) ;
}

struct BingwaComputation: public boost::noncopyable {
	typedef std::auto_ptr< BingwaComputation > UniquePtr ;
	static UniquePtr create( std::string const& name, std::vector< std::string > const& cohort_names, appcontext::OptionProcessor const& ) ;
	virtual ~BingwaComputation() {}
	typedef genfile::SNPIdentifyingData2 SNPIdentifyingData2 ;
	typedef boost::function< void ( std::string const& value_name, genfile::VariantEntry const& value ) > ResultCallback ;
	
	struct DataGetter: public boost::noncopyable {
		virtual ~DataGetter() {} ;
		virtual std::size_t get_number_of_cohorts() const = 0 ;
		virtual bool is_non_missing( std::size_t i ) const = 0 ;
		virtual void get_counts( std::size_t, Eigen::VectorXd* result ) const = 0 ;
		virtual void get_betas( std::size_t i, Eigen::VectorXd* result ) const = 0 ;
		virtual void get_ses( std::size_t i, Eigen::VectorXd* result  ) const = 0 ;
		virtual void get_covariance_upper_triangle( std::size_t i, Eigen::VectorXd* result  ) const = 0 ;
		virtual void get_pvalue( std::size_t i, double* result ) const = 0 ;
		virtual void get_info( std::size_t i, double* result ) const = 0 ;
		virtual void get_maf( std::size_t i, double* result ) const = 0 ;
		virtual void get_variable( std::string const& variable, std::size_t i, std::string* result ) const = 0 ;
	} ;

	typedef boost::function< bool ( DataGetter const&, int i ) > Filter ;
	
	virtual void set_effect_parameter_names( EffectParameterNamePack const& names ) = 0 ;
	virtual void get_variables( boost::function< void ( std::string ) > ) const = 0 ;
	virtual void operator()(
		SNPIdentifyingData2 const&,
		DataGetter const& data_getter,
		ResultCallback callback
	) = 0 ;
	virtual std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const = 0 ;
	virtual std::string get_spec() const = 0 ;
} ;

struct PerCohortValueReporter: public BingwaComputation {
public:
	typedef std::auto_ptr< PerCohortValueReporter > UniquePtr ;
	
public:
	PerCohortValueReporter( std::vector< std::string > const& cohort_names ):
		m_cohort_names( cohort_names )
	{}
	
	void add_variable( std::string const& variable ) {
		m_extra_variables.insert( variable ) ;
	}
	
	void set_effect_parameter_names( EffectParameterNamePack const& names ) {
		m_effect_parameter_names = names ;
	}
	
	void get_variables( boost::function< void ( std::string ) > callback ) const {
		using genfile::string_utils::to_string ;
		
		std::size_t const N = m_cohort_names.size() ;
		for( std::size_t i = 0; i < N; ++i ) {
			std::string prefix = m_cohort_names[ i ] + ":" ;
			callback( prefix + "A" ) ;
			callback( prefix + "B" ) ;
			callback( prefix + "AA" ) ;
			callback( prefix + "AB" ) ;
			callback( prefix + "BB" ) ;
			callback( prefix + "NULL" ) ;
			callback( prefix + "B_allele_frequency" ) ;
			callback( prefix + "maf" ) ;

			for( std::size_t i = 0; i < m_effect_parameter_names.size(); ++i ) {
//				callback( prefix + "beta_" + to_string( i+1 ) ) ;
//				callback( prefix + "se_" + to_string( i+1 ) ) ;
				callback( prefix + m_effect_parameter_names.parameter_name(i) ) ;
				callback( prefix + m_effect_parameter_names.se_name(i) ) ;
			}
			for( std::size_t i = 0; i < m_effect_parameter_names.size(); ++i ) {
				for( std::size_t j = i+1; j < m_effect_parameter_names.size(); ++j ) {
					callback( prefix + m_effect_parameter_names.covariance_name(i,j) ) ;
				}
			}
			
			callback( prefix + "pvalue" ) ;
			callback( prefix + "info" ) ;

			BOOST_FOREACH( std::string const& variable, m_extra_variables ) {
				callback( prefix + variable ) ;
			}
		}
	}
	
	void operator()(
		SNPIdentifyingData2 const&,
		DataGetter const& data_getter,
		ResultCallback callback
	) {
		std::size_t const N = m_cohort_names.size() ;
		for( std::size_t i = 0; i < N; ++i ) {
			if( data_getter.is_non_missing( i ) ) {
				Eigen::VectorXd betas ;
				Eigen::VectorXd ses ;
				Eigen::VectorXd covariance ;
				Eigen::VectorXd counts ;
				double pvalue ;
				double info ;
				double maf ;
				data_getter.get_betas( i, &betas ) ;
				data_getter.get_ses( i, &ses ) ;
				data_getter.get_covariance_upper_triangle( i, &covariance ) ;
				
				data_getter.get_counts( i, &counts ) ;
				data_getter.get_pvalue( i, &pvalue ) ;
				data_getter.get_info( i, &info ) ;
				data_getter.get_maf( i, &maf ) ;

				assert( counts.size() == 6 ) ;
				using genfile::string_utils::to_string ;
				std::string prefix = m_cohort_names[ i ] + ":" ;
				callback( prefix + "A", counts(0) ) ;
				callback( prefix + "B", counts(1) ) ;
				callback( prefix + "AA", counts(2) ) ;
				callback( prefix + "AB", counts(3) ) ;
				callback( prefix + "BB", counts(4) ) ;
				callback( prefix + "NULL", counts(5) ) ;

				double B_allele_count = 0 ;
				double total_allele_count = 0 ;
				if( counts(0) == counts(0) ) {
					B_allele_count += counts(1) ;
					total_allele_count += counts(0) + counts(1) ;
				}
				if( counts(2) == counts(2) ) {
					B_allele_count += counts(3) + 2 * counts(4) ;
					total_allele_count += 2.0 * ( counts(2) + counts(3) + counts(4) ) ;
				}
				callback( prefix + "B_allele_frequency", B_allele_count / total_allele_count ) ;
				callback( prefix + "maf", maf ) ;
				
				assert( betas.size() == ses.size() ) ;
				assert( covariance.size() == ( betas.size() - 1 ) * betas.size() / 2 ) ;
				for( int j = 0; j < betas.size(); ++j ) {
//					callback( prefix + "beta_" + to_string( j+1 ), betas(j) ) ;
					callback( prefix + m_effect_parameter_names.parameter_name(j), betas(j) ) ;
				}
				for( int j = 0; j < betas.size(); ++j ) {
					callback( prefix + m_effect_parameter_names.se_name(j), ses(j) ) ;
				}
				{
					int index = 0 ;
					for( int j = 0; j < betas.size(); ++j ) {
						for( int k = j+1; k < betas.size(); ++k, ++index ) {
							callback( prefix + m_effect_parameter_names.covariance_name(j,k), covariance( index ) ) ;
						}
					}
					assert( index == covariance.size() ) ;
				}
				callback( prefix + "pvalue", pvalue ) ;
				callback( prefix + "info", info ) ;
				{
					std::string value ;
					BOOST_FOREACH( std::string const& variable, m_extra_variables ) {
						data_getter.get_variable( variable, i, &value ) ;
						callback( prefix + variable, value ) ;
					}
				}
			}
		}
	}
	
	std::string get_spec() const {
		return "PerCohortValueReporter" ;
	}
	
	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
		return prefix + get_spec() ;
	}

	private:
		std::vector< std::string > const m_cohort_names ;
		std::set< std::string > m_extra_variables ;
		EffectParameterNamePack m_effect_parameter_names ;
} ;

namespace impl {
	bool basic_missingness_filter( BingwaComputation::DataGetter const& data_getter, int i ) {
		return data_getter.is_non_missing( i ) ;
	}
	
	bool info_maf_filter( BingwaComputation::DataGetter const& data_getter, int i, double const lower_info_threshhold, double const lower_maf_threshhold ) {
		bool result = data_getter.is_non_missing( i ) ;
		if( result ) {
			double info ;
			double maf ;
			data_getter.get_info( i, &info ) ;
			data_getter.get_maf( i, &maf ) ;
			if( (
				lower_info_threshhold == lower_info_threshhold && ( info != info || info < lower_info_threshhold )
			) || (
				lower_maf_threshhold == lower_maf_threshhold && ( maf != maf || maf < lower_maf_threshhold )
			) ) {
				result = false ;
			}
		}
		return result ;
	}
	
	bool get_betas_and_ses_one_per_study( BingwaComputation::DataGetter const& data_getter, BingwaComputation::Filter filter, Eigen::VectorXd& betas, Eigen::VectorXd& ses, Eigen::VectorXd& non_missingness ) {
		std::size_t N = data_getter.get_number_of_cohorts() ;
		betas.resize( N ) ;
		ses.resize( N ) ;
		non_missingness.resize( 2 * N ) ;
		for( std::size_t i = 0; i < N; ++i ) {
			non_missingness( i ) = filter( data_getter, i ) ? 1.0 : 0.0 ;

			if( non_missingness( i )) {
				Eigen::VectorXd data ;
				data_getter.get_betas( i, &data ) ;
				if( data.size() != 1 ) {
					return false ;
				}
				betas(i) = data(0) ;

				data_getter.get_ses( i, &data ) ;
				if( data.size() != 1 ) {
					return false ;
				}
				ses(i) = data(0) ;
				
				if( betas(i) != betas(i) || ses(i) != ses(i) ) {
					non_missingness(i) = 0 ;
				}
			}
		}
		return true ;
	}
	
	bool get_betas_and_covariance_per_study(
		BingwaComputation::DataGetter const& data_getter,
		BingwaComputation::Filter filter,
		Eigen::VectorXd& betas,
		Eigen::MatrixXd& covariance,
		Eigen::VectorXd& non_missingness,
		int const l
	) {
		std::size_t N = data_getter.get_number_of_cohorts() ;

		betas.setConstant( l*N, NA ) ;
		non_missingness.setConstant( l*N, 1 ) ;
		covariance.setZero( l*N, l*N ) ;

		for( std::size_t i = 0; i < N; ++i ) {
			non_missingness( i ) = filter( data_getter, i ) ? 1.0 : 0.0 ;

			if( non_missingness( i )) {
				Eigen::VectorXd this_betas, this_ses, this_covariance ;
				data_getter.get_betas( i, &this_betas ) ;
				data_getter.get_ses( i, &this_ses ) ;
				data_getter.get_covariance_upper_triangle( i, &this_covariance ) ;
				if( this_betas.size() != l || this_ses.size() != l || this_covariance.size() != ((l-1) * l / 2 ) ) {
					return false ;
				}
				// beta_1's for all pops go in one contiguous block.
				// Then beta_2's for all pops.
				// This matches the layout we use for the priors.
				int cov_i = 0 ;
				for( int j = 0; j < l; ++j ) {
					int const index = (j*N)+i; //(l*i)+j ;
					betas( index ) = this_betas(j) ;
					covariance( index, index ) = this_ses(j) * this_ses(j) ;
					if( this_betas(j) != this_betas(j) || this_ses(j) != this_ses(j) ) {
						non_missingness( index ) = 0 ;
					}
					for( int k = (j+1); k < l; ++k, ++cov_i ) {
						int const index2 = (k*N)+i; //(l*i)+j ;
						covariance( index, index2 ) = covariance( index2, index ) = this_covariance( cov_i ) ;
					}
				}
			}
		}
		return true ;
	}

#if 0
	bool get_fixed_effect_betas_and_covariance_from_study(
		BingwaComputation::DataGetter const& data_getter,
		BingwaComputation::Filter filter,
		Eigen::VectorXd& betas,
		Eigen::MatrixXd& covariance,
		Eigen::VectorXd& non_missingness,
		int const l
	) {
		std::size_t N = data_getter.get_number_of_cohorts() ;

		betas.setConstant( l, NA ) ;
		non_missingness.setConstant( l, 1 ) ;
		covariance.setZero( l, l ) ;

		// We use a temporary matrix A to accumulate the inverses of per-study covariances.
		Eigen::MatrixXd A = Eigen::MatrixXd::Zero( l, l ) ;
		
		for( std::size_t i = 0; i < N; ++i ) {
			non_missingness( i ) = filter( data_getter, i ) ? 1.0 : 0.0 ;

			if( non_missingness( i )) {
				Eigen::VectorXd this_betas, this_ses, this_covariance_upper_triangle ;
				data_getter.get_betas( i, &this_betas ) ;
				data_getter.get_ses( i, &this_ses ) ;
				data_getter.get_covariance_upper_triangle( i, &this_covariance_upper_triangle ) ;
				if( this_betas.size() != l || this_ses.size() != l || this_covariance_upper_triangle.size() != ((l-1) * l / 2 ) ) {
					return false ;
				}
				// Populate
				Eigen::MatrixXd const& this_covariance = Eigen::MatrixXd::Zero( l, l ) ;
				for( std::size_t i = 0; i < l; ++i ) {
					this_covariance.block( i, i, 1, l-i ) = this_covariance_upper_triangle.segment( (i*l) ) ;
				}

				// FIX ME
				FIX ME
			}
		}
		return true ;
	}
#endif
	std::string serialise( Eigen::VectorXd const& vector ) {
		std::ostringstream ostr ;
		for( int i = 0; i < vector.size(); ++i ) {
			if( i > 0 ) {
				ostr << "," ;
			}
			if( vector(i) == vector(i) ) {
				ostr << vector(i) ;
			} else {
				ostr << "NA" ;
			}
		}
		return ostr.str() ;
	}
}

/*
	Here is R code to do it:
	frequentist_meta_analysis <- function( betas, ses, side = NULL ) {
		if( class( betas ) == "numeric" ) {
			betas = matrix( betas, ncol = length( betas ))
			ses = matrix( ses, ncol = length( betas ))
		}
		v = ses^2 ;
		betas = betas / v ;
		denominator_terms = 1 / v;
		meta_beta = rowSums( betas ) / rowSums( denominator_terms )
		meta_se = sqrt( 1 / rowSums( denominator_terms ) )
		if( is.null( side ) ) {
			pvalue = 2 * pnorm( abs( meta_beta ), mean = 0, sd = meta_se, lower.tail = FALSE  ) ;
		} else {
			pvalue = 1 ;
			if( meta_beta * side > 0 ) {
				pvalue = pnorm( abs( meta_beta ), mean = 0, sd = meta_se, lower.tail = FALSE )
			}
		}
		return(
			list(
				meta_beta = meta_beta,
				meta_se = meta_se,
				pvalue = pvalue,
				side = side
			)
		) ;
	}
*/
struct FixedEffectFrequentistMetaAnalysis: public BingwaComputation {
	typedef std::auto_ptr< FixedEffectFrequentistMetaAnalysis > UniquePtr ;
	
	FixedEffectFrequentistMetaAnalysis():
		m_filter( &impl::basic_missingness_filter ),
		m_degrees_of_freedom( -1 )
	{}
	
	void set_filter( Filter filter ) {
		assert( filter ) ;
		m_filter = filter ;
	}
	
	void set_effect_parameter_names( EffectParameterNamePack const& names ) {
		m_degrees_of_freedom = names.size() ;
	}

	void get_variables( boost::function< void ( std::string ) > callback ) const {
		callback( "FixedEffectMetaAnalysis:meta_beta" ) ;
		callback( "FixedEffectMetaAnalysis:meta_se" ) ;
		callback( "FixedEffectMetaAnalysis:pvalue" ) ;
	}

	void operator()(
		SNPIdentifyingData2 const& snp,
		DataGetter const& data_getter,
		ResultCallback callback
	) {
		std::size_t const N = data_getter.get_number_of_cohorts() ;
		if( N == 0 ) {
			return ;
		}
		Eigen::VectorXd betas = Eigen::VectorXd::Constant( N, NA ) ;
		Eigen::VectorXd ses = Eigen::VectorXd::Constant( N, NA ) ;
		Eigen::VectorXd non_missingness = Eigen::VectorXd::Constant( N, NA ) ;
		if( !impl::get_betas_and_ses_one_per_study( data_getter, m_filter, betas, ses, non_missingness ) ) {
			callback( "FixedEffectMetaAnalysis:meta_beta", genfile::MissingValue() ) ;
			callback( "FixedEffectMetaAnalysis:meta_se", genfile::MissingValue() ) ;
			callback( "FixedEffectMetaAnalysis:pvalue", genfile::MissingValue() ) ;
			return ;
		}
		else {
			Eigen::VectorXd inverse_variances = ( ses.array() * ses.array() ).inverse() ;
			for( int i = 0; i < int(N); ++i ) {
				if( non_missingness( i ) == 0.0 ) {
					inverse_variances( i ) = 0 ;
					betas( i ) = 0 ;
					ses( i ) = 0 ;
				}
			}
			double const meta_beta = ( non_missingness.sum() == 0 ) ? NA : ( inverse_variances.array() * betas.array() ).sum() / inverse_variances.sum() ;
			double const meta_se = ( non_missingness.sum() == 0 ) ? NA : std::sqrt( 1.0 / inverse_variances.sum() ) ;

			callback( "FixedEffectMetaAnalysis:meta_beta", meta_beta ) ;
			callback( "FixedEffectMetaAnalysis:meta_se", meta_se ) ;

			//std::cerr << "SNP: " << snp << ": betas = " << betas << ", ses = " << ses << ".\n" ;

			if( meta_beta == meta_beta && meta_se == meta_se && meta_se > 0 && meta_se < std::numeric_limits< double >::infinity() ) {
				typedef boost::math::normal NormalDistribution ;
				NormalDistribution normal( 0, meta_se ) ;
				// P-value is the mass under both tails of the normal distribution larger than |meta_beta|
				double const pvalue = 2.0 * boost::math::cdf( boost::math::complement( normal, std::abs( meta_beta ) ) ) ;
				callback( "FixedEffectMetaAnalysis:pvalue", pvalue ) ;
			}
			
		}
	}

	std::string get_spec() const {
		return "FixedEffectFrequentistMetaAnalysis" ;
	}

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
		return prefix + get_spec() ;
	}
private:
	Filter m_filter ;
	int m_degrees_of_freedom ;
} ;


namespace impl {
	template< typename M >
	std::string format_matrix( M const& m ) {
		std::ostringstream s ;
		s << "matrix( nrow=" << m.rows() << ", ncol=" << m.cols() << ", data = c(" ;
		for( int i = 0; i < m.rows(); ++i ) {
			if( i > 0 ) {
				s << "," ;
			}
			for( int j = 0; j < m.cols(); ++j ) {
				if( j > 0 ) {
					s << "," ;
				}
				s << m(i,j) ;
			}
		}
		s << "))" ;
		return s.str() ;
	}

//
// Bayesian meta-analysis treating the prior as normal with mean 0 and variance matrix \Sigma
// and the likelihood as normal and given by
// the estimated betas and ses.
//
// Here is R code for the same computation:
/*
	bayesian_meta_analysis <- function( betas, ses, prior ) {
		if( length( ses ) > 1 ) {
			V = diag( ses^2 ) ;
		} else {
			V = matrix( data = ses^2, nrow = 1, ncol = 1 )
		}
		betas = matrix( betas, ncol = 1, nrow = length( betas )) ;
		constant = sqrt( det( V ) / det( V + prior ) ) ;
		exponent = 0.5 * t( betas ) %*% ( solve( V ) - solve( V + prior ) ) %*% betas
		print( V )
		print( betas )
		print( constant )
		print( exponent )
		return( constant * exp( exponent ))
	}
*/
	double compute_bayes_factor( Eigen::MatrixXd const& prior, Eigen::MatrixXd const& V, Eigen::VectorXd const& betas ) {
	#if DEBUG_BINGWA
		 std::cerr << std::resetiosflags( std::ios::floatfield ) ;
		 std::cerr << "prior = " << prior << ".\n" ;
		 std::cerr << "betas = " << betas.transpose() << ".\n" ;
		 std::cerr << "V = " << V << ".\n" ;
	#endif
		// I hope LDLT copes with noninvertible matrices.
		// Maybe it doesn't...but let's find out.
		Eigen::LDLT< Eigen::MatrixXd > Vsolver( V ) ;
		Eigen::LDLT< Eigen::MatrixXd > V_plus_prior_solver( V + prior ) ;
		Eigen::VectorXd exponent = betas.transpose() * ( Vsolver.solve( betas ) - V_plus_prior_solver.solve( betas ) ) ;
	
		assert( exponent.size() == 1 ) ;
	
		double const constant = std::sqrt( Vsolver.vectorD().prod() / V_plus_prior_solver.vectorD().prod() ) ;
		double const result = constant * std::exp( 0.5 * exponent(0) ) ;

	#if DEBUG_BINGWA
		std::cerr << "constant = " << constant << ".\n" ;
		std::cerr << "exponent= " << exponent.transpose() << ".\n" ;
		std::cerr << "BF = " << result << ".\n" ;
	#endif
		return result ;
	}


	Eigen::MatrixXd get_nonmissing_coefficient_selector(
		Eigen::VectorXd const& non_missingness
	) {
		int const number_of_included_effects = (
			( non_missingness.array() > 0 ).cast< double >()
		).sum() ;

		Eigen::MatrixXd result = Eigen::MatrixXd::Zero( number_of_included_effects, non_missingness.size() ) ;
	
		if( number_of_included_effects > 0 ) {
			int count = 0 ;
			for( int i = 0; i < non_missingness.size(); ++i ) {
				// std::cerr << "i=" << i << ", non_missingness(i) = " << non_missingness(i) << ", sigma( i, i ) = " << sigma(i,i) << ".\n" ;
				if( non_missingness(i) ) {
					result( count++, i ) = 1 ;
				}
			}
		
#if DEBUG_BINGWA	
			std::cerr << "impl::get_nonmissing_coefficient_selector()" << ": prior selector is:\n" << result << "\n" ;
#endif
		}
		return result ;
	}

}


struct ApproximateBayesianMetaAnalysis: public BingwaComputation {
	typedef std::auto_ptr< ApproximateBayesianMetaAnalysis > UniquePtr ;
	
	ApproximateBayesianMetaAnalysis(
		std::string const& name,
		Eigen::MatrixXd const& sigma
	):
		m_name( name ),
		m_prefix( name ),
		m_sigma( sigma ),
		m_filter( &impl::basic_missingness_filter )
	{
		assert( m_sigma.rows() == m_sigma.cols() ) ;
	}

	void set_filter( Filter filter ) {
		m_filter = filter ;
	}

	void set_effect_parameter_names( EffectParameterNamePack const& names ) {
		m_effect_parameter_names = names ;
	}

	void get_variables( boost::function< void ( std::string ) > callback ) const {
		callback( m_prefix + ":bf" ) ;
	}

	void operator()(
		SNPIdentifyingData2 const& snp,
		DataGetter const& data_getter,
		ResultCallback callback
	) {
		std::size_t const N = data_getter.get_number_of_cohorts() ;
		if( N == 0 ) {
			return ;
		}

		Eigen::VectorXd betas ;
		Eigen::MatrixXd covariance ;
		Eigen::VectorXd non_missingness ;
	
		if(
			!impl::get_betas_and_covariance_per_study( data_getter, m_filter, betas, covariance, non_missingness, m_effect_parameter_names.size() )
			|| non_missingness.sum() == 0
		) {
			return ;
		}
		else {
			Eigen::MatrixXd prior_selector = impl::get_nonmissing_coefficient_selector( non_missingness ) ;
			betas = prior_selector * betas ;
			covariance = prior_selector * covariance * prior_selector.transpose() ;
			Eigen::MatrixXd prior = prior_selector * m_sigma * prior_selector.transpose() ;
			
			double const bf = impl::compute_bayes_factor( prior, covariance, betas ) ;
			callback( m_prefix + ":bf", bf ) ;
		}
	}
		
	std::string get_spec() const {
		return "ApproximateBayesianMetaAnalysis( " + m_name + " ) with prior:\n" + genfile::string_utils::to_string( m_sigma ) ;
	}

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
		return prefix + get_spec() ;
	}
private:
	std::string const m_name ;
	std::string const m_prefix ;
	Eigen::MatrixXd const m_sigma ;
	Filter m_filter ;
	EffectParameterNamePack m_effect_parameter_names ;
} ;

struct ModelAveragingBayesFactorAnalysis: public BingwaComputation {
public:
	struct ModelSpec {
	public:
		ModelSpec( std::string const& name, Eigen::MatrixXd const& covariance, double weight ):
			m_name( name ),
			m_covariance( covariance ),
			m_weight( weight )
		{}
			
		ModelSpec( ModelSpec const& other ):
			m_name( other.m_name ),
			m_covariance( other.m_covariance ),
			m_weight( other.m_weight )
		{}

		ModelSpec& operator=( ModelSpec const& other ) {
			m_name = other.m_name ;
			m_covariance = other.m_covariance ;
			m_weight = other.m_weight ;
			return *this ;
		}
		
		std::string const& name() const { return m_name ; }
		Eigen::MatrixXd const& covariance() const { return m_covariance ; }
		double weight() const { return m_weight ; }
	private:
		std::string  m_name ;
		Eigen::MatrixXd m_covariance ;
		double m_weight ;
	} ;

public:
	ModelAveragingBayesFactorAnalysis() {}
	~ModelAveragingBayesFactorAnalysis() {}
	
	void add_model( std::string const& name, Eigen::MatrixXd const& sigma, double const& weight ) {
		assert( weight == weight ) ;
		m_models.push_back( ModelSpec( name, sigma, weight ) ) ;
	}
	
	void set_filter( Filter filter ) {
		m_filter = filter ;
	}

	void set_effect_parameter_names( EffectParameterNamePack const& names ) {
		m_effect_parameter_names = names ;
	}

	void get_variables( boost::function< void ( std::string ) > callback ) const {
		callback( m_prefix + ":bf" ) ;
	}

	void operator()(
		SNPIdentifyingData2 const& snp,
		DataGetter const& data_getter,
		ResultCallback callback
	) {
		std::size_t const N = data_getter.get_number_of_cohorts() ;
		if( N == 0 ) {
			return ;
		}

		Eigen::VectorXd betas ;
		Eigen::MatrixXd covariance ;
		Eigen::VectorXd non_missingness ;
	
		if(
			!impl::get_betas_and_covariance_per_study( data_getter, m_filter, betas, covariance, non_missingness, m_effect_parameter_names.size() )
			|| non_missingness.sum() == 0
		) {
			return ;
		}

		Eigen::MatrixXd prior_selector = impl::get_nonmissing_coefficient_selector( non_missingness ) ;
		betas = prior_selector * betas ;
		covariance = prior_selector * covariance ;
		
		std::vector< double > bfs( m_models.size(), NA ) ;
		std::vector< double > posteriors( m_models.size(), NA ) ;
		double total_weight ;
		double mean_bf = 0 ;
		double max_bf = -Infinity ;
		std::size_t max_bf_i = m_models.size() ;
		double max_posterior = -Infinity ;
		std::size_t max_posterior_i = m_models.size() ;

		for( std::size_t i = 0; i < m_models.size(); ++i ) {
			Eigen::MatrixXd const prior = prior_selector * m_models[i].covariance() * prior_selector.transpose() ;
			double const bf = impl::compute_bayes_factor( prior, covariance, betas ) ;
			if( bf == bf ) {
				double const weight = m_models[i].weight() ;
				total_weight += weight ;
				mean_bf += m_models[i].weight() * bf ;

				if( bf >= max_bf ) {
					max_bf = bf ;
					max_bf_i = i ;
				}

				double const posterior = weight * bf ;
				if( posterior >= max_posterior ) {
					max_posterior = posterior ;
					max_posterior_i = i ;
				}
			}
		}
		
		for( std::size_t i = 0; i < m_models.size(); ++i ) {
			if( bfs[i] == bfs[i] ) {
				callback( m_prefix + ":" + m_models[i].name() + ":bf", bfs[i] ) ;
			}
		}

		for( std::size_t i = 0; i < m_models.size(); ++i ) {
			if( bfs[i] == bfs[i] ) {
				callback( m_prefix + ":" + m_models[i].name() + ":posterior", posteriors[i] / total_weight ) ;
			}
		}
		
		callback( m_prefix + ":mean_bf", mean_bf / total_weight ) ;
		callback( m_prefix + ":max_bf", max_bf ) ;
		callback( m_prefix + ":max_bf_model", m_models[max_bf_i].name() ) ;
		callback( m_prefix + ":max_posterior", posteriors[ max_posterior_i ] / total_weight ) ;
		callback( m_prefix + ":max_posterior_model", m_models[ max_posterior_i ].name() ) ;
	}

	std::string get_spec() const {
		return "ModelAveragingBayesFactorAnalysis( " + m_name + " )" ;
	}

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
		return prefix + get_spec() ;
	}
private:
	std::string const m_name ;
	std::string const m_prefix ;
	Eigen::MatrixXd const m_sigma ;
	Filter m_filter ;
	EffectParameterNamePack m_effect_parameter_names ;
	

	std::vector< ModelSpec > m_models ;
} ;

BingwaComputation::UniquePtr BingwaComputation::create( std::string const& name, std::vector< std::string > const& cohort_names, appcontext::OptionProcessor const& options ) {
	BingwaComputation::UniquePtr result ;
	if( name == "FixedEffectFrequentistMetaAnalysis" ) {
		result.reset( new FixedEffectFrequentistMetaAnalysis() ) ;
	}
	else if( name == "PerCohortValueReporter" ) {
		PerCohortValueReporter::UniquePtr pcv( new PerCohortValueReporter( cohort_names ) ) ;
		if( options.check( "-extra-columns" )) {
			BOOST_FOREACH( std::string const& variable, options.get_values( "-extra-columns" )) {
				pcv->add_variable( variable ) ;
			}
		}
		result.reset( pcv.release() ) ;
	}
	else {
		throw genfile::BadArgumentError( "BingwaComputation::create()", "name=\"" + name + "\"" ) ;
	}
	return result ;
}

struct ValueAccumulator: public BingwaComputation {
public:
	typedef boost::shared_ptr< ValueAccumulator > SharedPtr ;
	typedef std::auto_ptr< ValueAccumulator > UniquePtr ;
public:
	ValueAccumulator( std::string const& name ):
		m_name( name ),
		m_accumulation( 0 ),
		m_accumulation_count( 0 )
	{}

	void set_weight( std::string const& model, double weight ) {
		std::cerr << "Setting weight for model " << model << " to " << weight << ".\n" ;
		m_weights[ model ] = weight ;
	}
	
	void set_effect_parameter_names( EffectParameterNamePack const& names ) {}

	void get_variables( boost::function< void ( std::string ) > callback ) const {
		callback( m_name ) ;
	}

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
		return prefix + get_spec()  ;
	}
	
	std::string get_spec() const {
		using genfile::string_utils::to_string ;
		std::string result =  "ValueAccumulator(" + m_name + ") with weights: " ;
		std::map< std::string, double >::const_iterator i = m_weights.begin() ;
		std::map< std::string, double >::const_iterator const end_i = m_weights.end() ;
		for( ; i != end_i; ++i ) {
			result += i->first + "=" + to_string( i->second ) + "; " ;
		}
		return result ;
	}
	
	void accumulate( std::string const& variable, genfile::VariantEntry const& value ) {
	//	std::cerr << "ValueAccumulator: looking at variable " << variable << " with value " << value << ".\n"  ;
		if( !value.is_missing() ) {
			std::map< std::string, double >::const_iterator where = m_weights.find( variable ) ;
			if( where != m_weights.end() ) {
	//			std::cerr << "ValueAccumulator: found!\n" ;
				m_accumulation += value.as< double >() * where->second ;
				++m_accumulation_count ;
			}
		}
	}

	void operator()(
		SNPIdentifyingData2 const& snp,
		DataGetter const&,
		ResultCallback callback
	) {
		if( m_accumulation_count > 0 ) {
			callback( m_name, m_accumulation / m_accumulation_count ) ;
		}
		m_accumulation = 0 ;
		m_accumulation_count = 0 ;
	}
	
	
	
private:
	std::string const m_name ;
	std::map< std::string, double > m_weights ;
	double m_accumulation ;
	double m_accumulation_count ;
} ;

namespace impl {
	void assign_counts( Eigen::MatrixXd* result, Eigen::MatrixXd::Index const col, genfile::VariantEntry const& AA, genfile::VariantEntry const& AB, genfile::VariantEntry const& BB, genfile::VariantEntry const& missing ) {
		assert( result ) ;
		assert( result->rows() >= 4 ) ;
		assert( result->cols() > col ) ;
		(*result)( 0, col ) = AA.is_missing() ? NA : AA.as< double >() ;
		(*result)( 1, col ) = AB.is_missing() ? NA : AB.as< double >() ;
		(*result)( 2, col ) = BB.is_missing() ? NA : BB.as< double >() ;
		(*result)( 3, col ) = missing.is_missing() ? NA : missing.as< double >() ;
	}
	
	void assign_vector_elt( Eigen::VectorXd* result, Eigen::VectorXd* non_missingness, Eigen::VectorXd::Index const elt, genfile::VariantEntry const& value ) {
		assert( result ) ;
		assert( non_missingness ) ;
		assert( result->size() == non_missingness->size() ) ;
		assert( result->size() > elt ) ;
		if( value.is_missing() ) {
			(*non_missingness)( elt ) = 0.0 ;
			(*result)( elt ) = NA ;
		}
		else {
			(*result)( elt ) = value.as< double >() ;
		}
	}
}

struct BingwaProcessor: public boost::noncopyable
{
public:
		typedef std::auto_ptr< BingwaProcessor > UniquePtr ;
public:
	static UniquePtr create( genfile::SNPIdentifyingData2::CompareFields const& compare_fields ) {
		return UniquePtr( new BingwaProcessor( compare_fields ) ) ;
	}
	
	BingwaProcessor( genfile::SNPIdentifyingData2::CompareFields const& compare_fields ):
		m_snps( compare_fields ),
		m_flip_alleles_if_necessary( false )
	{
	}
	
	void set_flip_alleles( void ) {
		m_flip_alleles_if_necessary = true ;
	}
	
	void add_cohort( std::string const& name, FrequentistGenomeWideAssociationResults::UniquePtr results ) {
		assert( results.get() ) ;
		m_cohort_names.push_back( name ) ;
		m_cohorts.push_back( results ) ;
	}
	
	std::size_t get_number_of_cohorts() const {
		return m_cohorts.size() ;
	}

	FrequentistGenomeWideAssociationResults const& get_cohort( std::size_t i ) const {
		return m_cohorts[i] ;
	}
	
	void summarise( appcontext::UIContext& ui_context ) {
		using genfile::string_utils::to_string ;

#if DEBUG_BINGWA
		ui_context.logger() << "================================================\n" ;
		ui_context.logger() << "Cohort summary:\n" ;
		for( std::size_t i = 0; i < m_cohorts.size(); ++i ) {
			ui_context.logger() << " - Cohort " << (i+1) << ":\n" ;
			ui_context.logger() << m_cohorts[i].get_summary( "   - " ) ;
			ui_context.logger() << "\n   - First few SNPs are:\n" ;
			for( std::size_t snp_i = 0; snp_i < std::min( std::size_t( 5 ), m_cohorts[i].get_number_of_SNPs() ); ++snp_i ) {
				double info ;
				double frequency ;
				m_cohorts[i].get_info( snp_i, &info ) ;
				ui_context.logger() << "     " << m_cohorts[i].get_SNP( snp_i ) << " (frequency = " << frequency << ", info = " << info << ")\n";
			}
		}
#endif
		
		ui_context.logger() << "\n================================================\n" ;
		ui_context.logger() << "SNP Categories:\n" ;
		ui_context.logger() << "  " ;
		for( std::size_t i = 0; i < m_cohorts.size(); ++i ) {
			ui_context.logger() << std::setw( 12 ) << ( "in scan " + to_string( i+1 )) ;
		}
		ui_context.logger() << "\n" ;
		for( CategoryCounts::const_iterator i = m_category_counts.begin(); i != m_category_counts.end(); ++i ) {
			ui_context.logger() << "  " ;
			for( std::size_t j = 0; j < m_cohorts.size(); ++j ) {	
				ui_context.logger() << std::setw(12) << ( i->first[j] ? '*' : ' ' ) ;
			}
			ui_context.logger() << ": " << i->second << "\n" ;
		}
		ui_context.logger() << " TOTAL: " << m_snps.size() << "\n" ;
		ui_context.logger() << "================================================\n" ;
	}
	
	void add_computation( std::string const& name, BingwaComputation::UniquePtr computation ) {
		m_computations.push_back( computation ) ;
		m_computations.back().set_effect_parameter_names( m_effect_parameter_names ) ;
	}

	void set_effect_parameter_names( EffectParameterNamePack const& names ) {
		m_effect_parameter_names = names ;
	}

	EffectParameterNamePack const get_effect_parameter_names() const {
		return m_effect_parameter_names ;
	}

	int const get_number_of_effect_parameters() const {
		return int( m_effect_parameter_names.size() ) ;
	}

	void get_variables( boost::function< void ( std::string const& ) > callback ) const {
		for( std::size_t i = 0; i < m_computations.size(); ++i ) {
			m_computations[i].get_variables( callback ) ;
		}
	}

	void setup( appcontext::UIContext& ui_context ) {
		unsafe_setup( ui_context ) ;
	}

	void process( appcontext::UIContext& ui_context ) {
		unsafe_process( ui_context ) ;
	}

	typedef boost::signals2::signal<
		void (
			genfile::SNPIdentifyingData2 const& snp,
			std::string const& value_name,
			genfile::VariantEntry const& value
		)
	> ResultSignal ;

	void send_results_to( ResultSignal::slot_type callback ) {
		m_result_signal.connect( callback ) ;
	}
	
private:
	std::vector< std::string > m_cohort_names ;
	boost::ptr_vector< FrequentistGenomeWideAssociationResults > m_cohorts ;
	EffectParameterNamePack m_effect_parameter_names ;
	boost::ptr_vector< BingwaComputation > m_computations ;
	struct SnpMatch {
	public:
		SnpMatch( std::size_t index_, bool flip_ ): index( index_ ), flip( flip_ ) {}
		SnpMatch(): index(0), flip( false ) {}
		SnpMatch( SnpMatch const& other ): index( other.index ), flip( other.flip ) {}
		SnpMatch& operator=( SnpMatch const& other ) { index = other.index ; flip = other.flip ; return *this ; }
	public:
		std::size_t index ;
		bool flip ;
	} ;
	typedef boost::optional< SnpMatch > OptionalSnpMatch ;
	typedef std::map<
		genfile::SNPIdentifyingData2,
		std::vector< OptionalSnpMatch >,
		genfile::SNPIdentifyingData2::CompareFields
	> SnpMap ;
	SnpMap m_snps ;
	typedef std::map< std::vector< bool >, std::size_t > CategoryCounts ;
	CategoryCounts m_category_counts ;
	bool m_flip_alleles_if_necessary ;
	
	ResultSignal m_result_signal ;

private:
	struct DataGetter: public BingwaComputation::DataGetter {
		DataGetter(
			boost::ptr_vector< FrequentistGenomeWideAssociationResults > const& cohorts,
			std::vector< OptionalSnpMatch >const& indices
		):
			m_cohorts( cohorts ),
			m_indices( indices )
		{}
		
		std::size_t get_number_of_cohorts() const { return m_cohorts.size() ; }
		
		void get_counts( std::size_t i, Eigen::VectorXd* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_counts( m_indices[i]->index, result ) ;
				if( m_indices[i]->flip ) {
					// Swap genotype counts for haploid:
					result->segment( 0, 2 ).reverseInPlace() ;
					// ...and diploid samples:
					result->segment( 2, 3 ).reverseInPlace() ;
				}
			}
		}
		void get_betas( std::size_t i, Eigen::VectorXd* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_betas( m_indices[i]->index, result ) ;
				if( m_indices[i]->flip ) {
					(*result) *= -1 ;
				}
			}
		}
		void get_ses( std::size_t i, Eigen::VectorXd* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_ses( m_indices[i]->index, result ) ;
			}
		}
		void get_covariance_upper_triangle( std::size_t i, Eigen::VectorXd* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_covariance_upper_triangle( m_indices[i]->index, result ) ;
			}
		}
		void get_pvalue( std::size_t i, double* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_pvalue( m_indices[i]->index, result ) ;
			}
		}
		void get_info( std::size_t i, double* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_info( m_indices[i]->index, result ) ;
			}
		}
		void get_maf( std::size_t i, double* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_maf( m_indices[i]->index, result ) ;
			}
		}
		bool is_non_missing( std::size_t i ) const {
			return( m_indices[i] ) ;
		}
		void get_variable( std::string const& variable, std::size_t i, std::string* result ) const {
			m_cohorts[i].get_variable( m_indices[i]->index, variable, result ) ;
		}

		private:
			boost::ptr_vector< FrequentistGenomeWideAssociationResults > const& m_cohorts ;
			std::vector< OptionalSnpMatch > const& m_indices ;
	} ;

private:
	void unsafe_setup( appcontext::UIContext& ui_context ) {
		link_data( ui_context ) ;
		categorise_by_missingness() ;
	}
	
	void add_SNP_callback( std::size_t cohort_i, std::size_t snp_i, genfile::SNPIdentifyingData2 const& snp ) {
		// Find the SNP that matches the given one (if it exists)
		std::pair< SnpMap::iterator, SnpMap::iterator > range = m_snps.equal_range( snp ) ;
		SnpMatch snp_match( snp_i, false ) ;

		if( range.second == range.first && m_flip_alleles_if_necessary ) {
			genfile::SNPIdentifyingData2 swapped_snp = snp ;
			swapped_snp.swap_alleles() ;
			range = m_snps.equal_range( swapped_snp ) ;
			snp_match.flip = ( range.second != range.first ) ;
		}
		
		if( range.second == range.first ) {
			// no match, so add this SNP.
			std::vector< OptionalSnpMatch >& snp_matches = m_snps[ snp ] ;
			snp_matches.resize( m_cohorts.size() ) ;
			snp_matches[ cohort_i ] = snp_match ;
		}
		else {
			// There is a match.  We make sure to record all the IDs
			// presented for this variant, but take the first rsid seen
			// as the rsid.
			// First save the currently-stored data and erase the current record.
			genfile::SNPIdentifyingData2 stored_snp = range.first->first ;
			std::vector< OptionalSnpMatch > snp_matches = range.first->second ;
			m_snps.erase( range.first ) ;
			
			stored_snp.add_identifier( snp.get_rsid() ) ;
			snp.get_alternative_identifiers( boost::bind( &genfile::SNPIdentifyingData2::add_identifier, &stored_snp, _1 ) ) ;

			snp_matches[ cohort_i ] = snp_match ;
			m_snps[ stored_snp ] = snp_matches ;
		}
	}

	void link_data( appcontext::UIContext& ui_context ) {
		appcontext::UIContext::ProgressContext progress_context = ui_context.get_progress_context( "Linking SNPs" ) ;
		progress_context( 0, m_cohorts.size() ) ;
		for( std::size_t cohort_i = 0; cohort_i < m_cohorts.size(); ++cohort_i ) {
			m_cohorts[ cohort_i].get_SNPs(
				boost::bind(
					&BingwaProcessor::add_SNP_callback,
					this,
					cohort_i,
					_1,
					_2
				)
			) ;
			progress_context( cohort_i+1, m_cohorts.size() ) ;
		}
	}
	
	void categorise_by_missingness() {
		SnpMap::const_iterator snp_i = m_snps.begin(), end_i = m_snps.end() ;
		for( ; snp_i != end_i; ++snp_i ) {
			std::vector< OptionalSnpMatch > const& indices = snp_i->second ;
			std::vector< bool > indicator( indices.size(), false ) ;
			for( std::size_t i = 0; i < indices.size(); ++i ) {
				indicator[i] = bool( indices[i] ) ;
			}
			++m_category_counts[ indicator ] ;
		}
	}

	void unsafe_process( appcontext::UIContext& ui_context ) {
		Eigen::MatrixXd cohort_counts( 4, m_cohorts.size() ) ;
		Eigen::MatrixXd cohort_betas( 2, m_cohorts.size() ) ;
		Eigen::MatrixXd cohort_ses( 2, m_cohorts.size() ) ;
		Eigen::VectorXd non_missingness( m_cohorts.size() ) ;

		{
			appcontext::UIContext::ProgressContext progress_context = ui_context.get_progress_context( "Storing meta-analysis results" ) ;
			SnpMap::const_iterator snp_i = m_snps.begin() ;
			SnpMap::const_iterator const end_i = m_snps.end() ;
			for( std::size_t snp_index = 0; snp_i != end_i; ++snp_i, ++snp_index ) {
				std::vector< OptionalSnpMatch > const& indices = snp_i->second ;
				DataGetter data_getter( m_cohorts, indices ) ;
				for( std::size_t i = 0; i < m_computations.size(); ++i ) {
					m_computations[i](
						snp_i->first,
						data_getter,
						boost::bind(
							boost::ref( m_result_signal ),
							snp_i->first,
							_1, _2
						)
					) ;
				}

				progress_context( snp_index + 1, m_snps.size() ) ;
			}
		}
	}
} ;


struct BingwaApplication: public appcontext::ApplicationContext {
public:
	typedef std::map< std::string, std::vector< std::string > > GroupDefinition ;
	typedef std::map< std::string, std::vector< std::string > > ValueListSet ;

public:
	BingwaApplication( int argc, char **argv ):
		appcontext::ApplicationContext(
			globals::program_name,
			globals::program_version,
			std::auto_ptr< appcontext::OptionProcessor >( new BingwaOptions ),
			argc,
			argv,
			"-log"
		),
		m_processor( BingwaProcessor::create( genfile::SNPIdentifyingData2::CompareFields( options().get< std::string > ( "-snp-match-fields" )) ) )
	{
		if( options().check( "-flip-alleles" )) {
			m_processor->set_flip_alleles() ;
		}
	}
	
	void run() {
		try {
			unsafe_run() ;
		}
		catch( genfile::InputError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): " << e.format_message() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( genfile::FileNotFoundError const& e ) {
			get_ui_context().logger() << "!! Error: No file matching \"" << e.filespec() << "\" could be found.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( genfile::OperationFailedError const& e ) {
			get_ui_context().logger() << "!! Error: operation failed: " << e.get_message() << "\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( db::Error const& e ) {
			get_ui_context().logger() << "!! Database error (" << e.what() << "): " << e.description() << ".\n" ;
			get_ui_context().logger() << "!! Error code is " << e.error_code() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
	}

	void unsafe_run() {
		load_data() ;

		m_processor->setup( get_ui_context() ) ;
		m_processor->summarise( get_ui_context() ) ;
		
		qcdb::Storage::SharedPtr storage ;
		if( options().check( "-flat-file" )) {
			storage = qcdb::FlatFileOutputter::create_shared(
				options().get< std::string >( "-o" ),
				options().get< std::string >( "-analysis-name" ),
				options().get_values_as_map()
			) ;
		}
		else {
			qcdb::FlatTableDBOutputter::SharedPtr table_storage = qcdb::FlatTableDBOutputter::create_shared(
				options().get< std::string >( "-o" ),
				options().get< std::string >( "-analysis-name" ),
				options().get< std::string >( "-analysis-chunk" ),
				options().get_values_as_map()
			) ;

			if( options().check( "-table-name" ) ) {
				table_storage->set_table_name( options().get< std::string >( "-table-name" )) ;
			}
			storage = table_storage ;
		}

		m_processor->send_results_to(
			boost::bind(
				&qcdb::Storage::store_per_variant_data,
				storage,
				_1, _2, _3
			)
		) ;

		std::vector< std::string > cohort_names ;
		if( options().check( "-cohort-names" )) {
			cohort_names = options().get_values< std::string >( "-cohort-names" ) ;
			if( cohort_names.size() != options().get_values< std::string >( "-data" ).size() ) {
				throw genfile::BadArgumentError(
					"BingwaApplication::unsafe_run()", "-cohort-names=\"" + genfile::string_utils::join( cohort_names, " " ) + "\"",
					"Number of cohort names does not match number of data files"
				) ;
			}
		}
		else {
			cohort_names.resize( options().get_values< std::string >( "-data" ).size() ) ;
			for( std::size_t i = 0; i < cohort_names.size(); ++i ) {
				cohort_names[i] = "cohort " + genfile::string_utils::to_string( i + 1 ) ;
			}
		}
	
		if( options().check( "-data" ) && options().get_values< std::string >( "-data" ).size() > 0 ) {
			m_processor->add_computation(
				"PerCohortValueReporter",
				BingwaComputation::create( "PerCohortValueReporter", cohort_names, options() )
			) ;
		
			if( options().check( "-no-meta-analysis" )) {
				// No more computations.
			}
			else {
				BingwaComputation::Filter filter = get_filter( options() ) ;

				{
					FixedEffectFrequentistMetaAnalysis::UniquePtr computation(
						new FixedEffectFrequentistMetaAnalysis()
					) ;
					computation->set_filter( filter ) ;
					m_processor->add_computation(
						"FixedEffectFrequentistMetaAnalysis",
						BingwaComputation::UniquePtr( computation.release() )
					) ;
				}
				if( options().check( "-simple-prior" ) || options().check( "-complex-prior" ) || options().check( "-group-specific-prior" )) {
					std::map< std::string, Eigen::MatrixXd > const priors = get_priors( options(), cohort_names ) ;
					assert( priors.size() > 0 ) ;
					{
						std::map< std::string, Eigen::MatrixXd >::const_iterator i = priors.begin() ;
						std::map< std::string, Eigen::MatrixXd >::const_iterator const end_i = priors.end() ;
						for( ; i != end_i; ++i ) {
							ApproximateBayesianMetaAnalysis::UniquePtr computation(
								new ApproximateBayesianMetaAnalysis(
									"ApproximateBayesianMetaAnalysis:" + i->first,
									i->second
								)
							) ;
						
							computation->set_filter( filter ) ;

							m_processor->add_computation(
								"ApproximateBayesianMetaAnalysis",
								BingwaComputation::UniquePtr( computation.release() )
							) ;
						}
					}
					
					// Now set up an average bayes factor
					// Currently this must be added after the computations above in order
					// that it is fed data (via the accumulate() callback) before the mean is
					// computed.
					{
						ValueAccumulator::UniquePtr mean_bf(
							new ValueAccumulator( "ApproximateBayesianMetaAnalysis:mean_bf" )
						) ;

						std::map< std::string, Eigen::MatrixXd >::const_iterator i = priors.begin() ;
						std::map< std::string, Eigen::MatrixXd >::const_iterator const end_i = priors.end() ;
						for( ; i != end_i; ++i ) {
							mean_bf->set_weight( "ApproximateBayesianMetaAnalysis:" + i->first + ":bf", 1.0 ) ;
						}

						m_processor->send_results_to(
							boost::bind(
								&ValueAccumulator::accumulate,
								mean_bf.get(),
								_2,
								_3
							)
						) ;
						
						m_processor->add_computation(
							"ApproximateBayesianMetaAnalysis",
							BingwaComputation::UniquePtr( mean_bf.release() )
						) ;
					}
				
					summarise_priors( priors, cohort_names ) ;
				}
			}
		}
		
		m_processor->get_variables(
			boost::bind(
				&qcdb::Storage::add_variable,
				storage,
				_1
			)
		) ;
		
		m_processor->process( get_ui_context() ) ;
		
		get_ui_context().logger() << "Finalising storage...\n" ;

		long finalise_options = ( options().check( "-noindex" ) ? 0 : qcdb::eCreateIndices ) ;
		storage->finalise( finalise_options ) ;
	}
	
	BingwaComputation::Filter get_filter( appcontext::OptionProcessor const& options ) const {
		BingwaComputation::Filter filter( &impl::basic_missingness_filter ) ;
		if( options.check( "-min-info" ) || options.check( "-min-maf" )) {
			double lower_info_threshhold = NA ;
			double lower_maf_threshhold = NA ;
			if( options.check( "-min-info" ) ) {
				lower_info_threshhold = options.get< double > ( "-min-info" ) ;
			}
			if( options.check( "-min-maf" ) ) {
				lower_maf_threshhold = options.get< double > ( "-min-maf" ) ;
			}
			filter = boost::bind(
				&impl::info_maf_filter, _1, _2,
				lower_info_threshhold, lower_maf_threshhold
			) ;
		}
		return filter ;
	}

	typedef boost::tuple< std::vector< std::string >, std::vector< std::string >, std::vector< std::string > > CovarianceSpec ;
	enum { eBetweenPopulationCorrelation = 0, eSD = 1, eCorrelation = 2 } ;

	void get_simple_priors(
		int const number_of_cohorts,
		std::vector< std::string > const& model_specs,
		boost::optional< std::vector< std::string > > model_names,
		std::map< std::string, Eigen::MatrixXd >* result
	) const {
		using namespace genfile::string_utils ;
		if( model_names && model_names.get().size() != model_specs.size() ) {
			throw genfile::BadArgumentError(
				"BingwaProcessor::get_simple_priors()",
				"model_names=\"" + join( model_names.get(), " " ) + "\"",
				"Wrong number of simple prior names specified (" + to_string( model_names.get().size() ) + ", should be " + to_string( model_specs.size() ) + ")"
			) ;
		}
		
		for( std::size_t i = 0; i < model_specs.size(); ++i ) {
			std::string const model_name = ( model_names ? model_names.get()[ i ] + "/" : "" ) ;
			std::vector< CovarianceSpec > const covariance_specs = parse_covariance_spec( model_specs[i] ) ;
			for( std::size_t j = 0; j < covariance_specs.size(); ++j ) {
				CovarianceSpec const& covariance_spec = covariance_specs[j] ;
				if( covariance_spec.get< eBetweenPopulationCorrelation >().size() > 0 ) {
					get_uniform_population_prior(
						number_of_cohorts,
						model_specs[ i ],
						covariance_spec,
						model_name,
						result
					) ;
				} else {
					get_all_population_prior(
						number_of_cohorts,
						model_specs[ i ],
						covariance_spec,
						model_name,
						result
					) ;
				}
			}
		}
	}

	void get_uniform_population_prior(
		int const number_of_cohorts,
		std::string const& model_spec,
		CovarianceSpec const& covariance_spec,
		std::string model_name,
		std::map< std::string, Eigen::MatrixXd >* result
	) const {
		using namespace genfile::string_utils ;
		
#if DEBUG_BINGWA
		std::cerr << "rho = " ;
		for( std::size_t j = 0; j < covariance_spec.get<eBetweenPopulationCorrelation>().size(); ++j ) {
			std::cerr << "\"" + covariance_spec.get<eBetweenPopulationCorrelation>()[j] << "\" " ;
		}
		std::cerr << "sd = " ;
		for( std::size_t j = 0; j < covariance_spec.get<eSD>().size(); ++j ) {
			std::cerr << "\"" + covariance_spec.get<eSD>()[j] << "\" " ;
		}
		std::cerr << "cor = " ;
		for( std::size_t j = 0; j < covariance_spec.get<eCorrelation>().size(); ++j ) {
			std::cerr << "\"" + covariance_spec.get<eCorrelation>()[j] << "\" " ;
		}
		std::cerr << std::endl ;
#endif

		try {
			if( covariance_spec.get<eBetweenPopulationCorrelation>().size() != 1 ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_simple_priors()",
					"model_spec=\"" + model_spec + "\"",
					"In model spec \"" + model_spec + "\", number of correlations specified ("
					+ to_string( covariance_spec.get<eBetweenPopulationCorrelation>().size() )
					+ ") is not 1"
				) ;
			}
			if( int( covariance_spec.get<eSD>().size() ) != m_processor->get_number_of_effect_parameters() ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_simple_priors()",
					"model_spec=\"" + model_spec + "\"",
					"In model spec \"" + model_spec + "\", number of sds specified ("
					+ to_string( covariance_spec.get<eSD>().size() )
					+ ") does not match number of effect size parameters (" + to_string( m_processor->get_number_of_effect_parameters() ) + ")"
				) ;
			}

			model_name += "rho=" + join( covariance_spec.get<eBetweenPopulationCorrelation>(), "," ) ;

			{
				int const N = m_processor->get_number_of_cohorts() ;
				int const D = m_processor->get_number_of_effect_parameters() ;
				Eigen::MatrixXd correlation = Eigen::MatrixXd::Zero( N * D, N * D ) ;

				// We start by computing the correlation matrix.
				// Prior is layed out in order of parameter (not cohort).
				// Thus the first NxN block is for the first effect size parameter.
				// We start with the diagonal blocks.  These do not involve the between-parameter correlation.
				for( int block_i = 0; block_i < D; ++block_i ) {
					correlation.block( block_i * N, block_i * N, N, N ) = get_prior_matrix(
						number_of_cohorts,
						to_repr< double >( covariance_spec.get<eBetweenPopulationCorrelation>()[0] ),
						1
					) ;
				}
				// ...and move on to the off-diagonal blocks which involve between-cohort and between-parameter correlation.
				// We compute this using the diagonal block
				int c = 0 ;
				for( int block_i = 0; block_i < D; ++block_i ) {
					if( covariance_spec.get<eCorrelation>()[ c++ ] != "1" ) {
						throw genfile::BadArgumentError(
							"BingwaProcessor::get_simple_priors()",
							"model_spec=\"" + model_spec + "\"",
							"Between-parameter correlation matrix must have 1's on diagonal"
						) ;
					}
					for( int block_j = block_i + 1; block_j < D; ++block_j, ++c ) {
						correlation.block( block_i * N, block_j * N, N, N ).array()
							= get_prior_matrix(
								number_of_cohorts,
								to_repr< double >( covariance_spec.get<eBetweenPopulationCorrelation>()[0] ),
								1
							) * to_repr< double >( covariance_spec.get<eCorrelation>()[c] ) ;
						;
						correlation.block( block_j * N, block_i * N, N, N ) = correlation.block( block_i * N, block_j * N, N, N ) ;
					}
				}

				Eigen::VectorXd sd_matrix( N * D ) ;
				for( int block_i = 0; block_i < D; ++block_i ) {
					sd_matrix.segment( block_i * N, N ).setConstant( to_repr< double >( covariance_spec.get<eSD>()[block_i] ) ) ;
				}
				Eigen::MatrixXd const covariance = sd_matrix.asDiagonal() * correlation * sd_matrix.asDiagonal() ;

				(*result)[ model_name + "/sd=" + join( covariance_spec.get<eSD>(), "," ) + "/cor=" + join( covariance_spec.get<eCorrelation>(), "," ) ] = covariance ;
				
#if DEBUG_BINGWA
				std::cerr << "For model " << ( model_name + "/sd=" + join( covariance_spec.get<eSD>(), " " ) ) << ":\n"
					<< "correlation =\n" << correlation << "\n"
					<< "covariance = \n" << covariance << "\n" ;
#endif
			}
		}
		catch( StringConversionError const& e ) {
			throw genfile::BadArgumentError(
				"BingwaProcessor::get_simple_priors()",
				"model_spec=\"" + model_spec + "\"",
				"Expected a numerical value."
			) ;
		}
	}

	void get_all_population_prior(
		int const number_of_cohorts,
		std::string const& model_spec,
		CovarianceSpec const& covariance_spec,
		std::string model_name,
		std::map< std::string, Eigen::MatrixXd >* result
	) const {
		using namespace genfile::string_utils ;

		int const N = m_processor->get_number_of_cohorts() ;
		int const D = m_processor->get_number_of_effect_parameters() ;
		int const dimension = N*D ;

#if DEBUG_BINGWA
		std::cerr << "rho = " ;
		for( std::size_t j = 0; j < covariance_spec.get<eBetweenPopulationCorrelation>().size(); ++j ) {
			std::cerr << "\"" + covariance_spec.get<eBetweenPopulationCorrelation>()[j] << "\" " ;
		}
		std::cerr << "sd = " ;
		for( std::size_t j = 0; j < covariance_spec.get<eSD>().size(); ++j ) {
			std::cerr << "\"" + covariance_spec.get<eSD>()[j] << "\" " ;
		}
		std::cerr << "cor = " ;
		for( std::size_t j = 0; j < covariance_spec.get<eCorrelation>().size(); ++j ) {
			std::cerr << "\"" + covariance_spec.get<eCorrelation>()[j] << "\" " ;
		}
		std::cerr << std::endl ;
#endif

		try {
			if( covariance_spec.get<eBetweenPopulationCorrelation>().size() != 0 ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_full_prior()",
					"model_spec=\"" + model_spec + "\"",
					"In model spec \"" + model_spec + "\", number of between-cohort correlations specified ("
					+ to_string( covariance_spec.get<eBetweenPopulationCorrelation>().size() )
					+ ") is not 0"
				) ;
			}
			if( int( covariance_spec.get<eSD>().size() ) != dimension ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_full_prior()",
					"model_spec=\"" + model_spec + "\"",
					"In model spec \"" + model_spec + "\", number of sds specified ("
					+ to_string( covariance_spec.get<eSD>().size() )
					+ ") does not match number of effect size parameters (" + to_string( dimension ) + ")"
				) ;
			}

			std::size_t const numberOfCorrelations = ( dimension * ( dimension + 1 ) / 2 ) ;
			if( int( covariance_spec.get<eCorrelation>().size() ) != numberOfCorrelations ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_full_prior()",
					"model_spec=\"" + model_spec + "\"",
					"In model spec \"" + model_spec + "\", number of correlations specified ("
					+ to_string( covariance_spec.get<eCorrelation>().size() )
					+ ") does not match number of effect size parameters (" + to_string( numberOfCorrelations ) + ")"
				) ;
			}

			{
				Eigen::MatrixXd correlation = Eigen::MatrixXd::Zero( dimension, dimension ) ;
				// We start by parsing the correlation matrix into the matrix.
				for( int c = 0, i = 0; i < correlation.rows(); ++i ) {
					for( int j = i; j < correlation.cols(); ++j, ++c ) {
						correlation(i,j) = to_repr< double >( covariance_spec.get<eCorrelation>()[c] ) ;
						if( j > i ) {
							correlation(j,i) = correlation(i,j) ;
						}
					}
				}

				Eigen::VectorXd sd_vector( dimension ) ;
				for( int i = 0; i < dimension; ++i ) {
					sd_vector(i) = to_repr< double >( covariance_spec.get<eSD>()[ i ] ) ;
				}

				Eigen::MatrixXd const covariance = sd_vector.asDiagonal() * correlation * sd_vector.asDiagonal() ;

				(*result)[ model_name + "sd=" + join( covariance_spec.get<eSD>(), "," ) + "/cor=" + join( covariance_spec.get<eCorrelation>(), "," ) ] = covariance ;
				
#if DEBUG_BINGWA
				std::cerr << "For model " << ( model_name + "/sd=" + join( covariance_spec.get<eSD>(), " " ) ) << ":\n"
					<< "correlation =\n" << correlation << "\n"
					<< "covariance = \n" << covariance << "\n" ;
#endif
			}
		}
		catch( StringConversionError const& e ) {
			throw genfile::BadArgumentError(
				"BingwaProcessor::get_full_prior()",
				"model_spec=\"" + model_spec + "\"",
				"Expected a numerical value."
			) ;
		}
	}
	
	void get_complex_priors(
		int const number_of_cohorts,
		std::vector< std::string > const& model_specs,
		std::vector< std::string > const& model_names,
		std::map< std::string, Eigen::MatrixXd >* result
	) {
		using namespace genfile::string_utils ;
		
		if( model_specs.size() != model_names.size() ) {
			throw genfile::BadArgumentError(
				"BingwaProcessor::get_complex_priors()",
				"model_names=\"" + join( model_names, " " ) + "\"",
				"Wrong number of complex prior names specified (" + to_string( model_names.size() ) + ", should be " + to_string( model_specs.size() ) + ")"
			) ;
		}

		for( std::size_t i = 0; i < model_specs.size(); ++i ) {
			std::vector< std::string > bits = split_and_strip( model_specs[i], "*", " \t\n" ) ;
			std::string sd = "1" ;
			if( bits.size() > 2 ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_complex_priors()",
					"model_spec #" + to_string( i + 1 ) + ": \"" + model_specs[i] + "\"",
					"Malformatted model spec, should be of the form [sd]^2 * [matrix]."
				) ;
			} else if( bits.size() == 2 ) {
				if( bits[0].size() < 3 || bits[0].substr( bits[0].size() - 2, 2 ) != "^2" ) {
					throw genfile::BadArgumentError(
						"BingwaProcessor::get_complex_priors()",
						"model_spec #" + to_string( i + 1 ) + ": \"" + model_specs[i] + "\"",
						"Malformatted model spec, should be of the form [sd]^2 * [matrix]."
					) ;
				}
				sd = bits[0].substr( 0, bits[0].size() - 2 ) ;
			}
			try {
				int const N = m_processor->get_number_of_cohorts() ;
				int const D = m_processor->get_number_of_effect_parameters() ;
				std::vector< std::string > sds = expand_sd( sd, m_value_sets[ "sd" ] ) ;
				std::cerr << ">>>> sds.size() = " << sds.size() << ".\n" ;
				for( std::size_t sd_i = 0; sd_i < sds.size(); ++sd_i ) {
					std::cerr << ">>>> Adding complex prior for sd = " << sds[sd_i] << ".\n" ;
					(*result)[ model_names[i] + "/sd=" + sds[sd_i] ] = get_prior_matrix( 
						N*D,
						bits.back(),
						to_repr< double >( sds[sd_i] )
					) ;
				}
			} catch( genfile::string_utils::StringConversionError const& e ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_complex_priors()",
					"model_spec #" + to_string( i + 1 ) + ": \"" + model_specs[i] + "\"",
					"Malformatted model spec, should be of the form [rho]^2 * [matrix]."
				) ;
			}
		}
	}

	std::vector< std::string > get_group_members( std::string const& group_name, GroupDefinition const& groups ) const {
		GroupDefinition::const_iterator where = groups.find( group_name ) ;
		if( where == groups.end() ) {
			throw genfile::BadArgumentError(
				"BingwaProcessor::get_group_members()",
				"group_name=\"" + group_name + "\"",
				"Group name \"" + group_name + "\" is not recognised, please define cohort groups using -groups."
			) ;
		}
		return where->second ;
	}
	
	std::vector< int > get_cohort_indices( std::vector< std::string > cohorts, std::vector< std::string > const& all_cohorts ) const {
		std::vector< int > result( cohorts.size() ) ;
		for( std::size_t k = 0; k < cohorts.size(); ++k ) {
			std::vector< std::string >::const_iterator const cohort_i = std::find( all_cohorts.begin(), all_cohorts.end(), cohorts[k] ) ;
			if( cohort_i == all_cohorts.end() ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_cohort_indices()",
					"cohort=\"" + cohorts[k] + "\"",
					"\"" + cohorts[k] + "\" is not a cohort name."
				) ;
			}
			result[k] = int( cohort_i - all_cohorts.begin() ) ;
		}
		return result ;
	}
	
	void get_group_priors(
		int const number_of_cohorts,
		GroupDefinition const& groups,
		std::vector< std::string > const& model_specs,
		boost::optional< std::vector< std::string > > model_names,
		std::vector< std::string > const& cohort_names,
		std::map< std::string, Eigen::MatrixXd >* result
	) {
		using namespace genfile::string_utils ;

		if( model_names && model_specs.size() != model_names.get().size() ) {
			throw genfile::BadArgumentError(
				"BingwaProcessor::get_group_priors()",
				"model_names=\"" + join( model_names.get(), " " ) + "\"",
				"Wrong number of group prior names specified (" + to_string( model_names.get().size() ) + ", should be " + to_string( model_specs.size() ) + ")"
			) ;
		}
		
		//
		// format is group1:rho=r1/sd=s1,group2:rho=r2/sd=s2,...:rho=0.2
		//
		for( std::size_t i = 0; i < model_specs.size(); ++i ) {
			std::vector< std::string > bits ;
			{
				std::size_t pos = model_specs[i].rfind( ":" ) ;
				if( pos == std::string::npos ) {
					throw genfile::BadArgumentError(
						"BingwaProcessor::get_group_priors()",
						"model_spec=\"" + model_specs[i] + "\"",
						"Model spec \"" + model_specs[i] + "\" is malformed."
					) ;
				}
				bits.push_back( model_specs[i].substr( 0, pos ) ) ;
				bits.push_back( model_specs[i].substr( pos+1, model_specs[i].size() ) ) ;
			}
	
			std::vector< std::string > const between_cohort_rho_string = parse_rho( bits[1] ) ;
			if( between_cohort_rho_string.size() != m_processor->get_number_of_effect_parameters() ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_group_priors()",
					"model_spec=\"" + model_specs[i] + "\"",
					"In model spec \"" + model_specs[i] + "\", number of correlations specified ("
					+ to_string( between_cohort_rho_string.size() )
					+ ") does not match number of effect size parameters (" + to_string( m_processor->get_number_of_effect_parameters() ) + ")"
				) ;
			}
			assert( between_cohort_rho_string.size() == 1 ) ;
			double between_cohort_rho = to_repr< double >( between_cohort_rho_string[0] ) ;
			
			std::vector< std::string > const group_specs = split_and_strip( bits[0], "," ) ;
			if( group_specs.size() == 0 ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_group_priors()",
					"model_spec=\"" + model_specs[i] + "\"",
					"In model spec \"" + model_specs[i] + "\", no groups are specified."
				) ;
			}
			
			Eigen::MatrixXd prior( number_of_cohorts, number_of_cohorts ) ;
			prior.setZero() ;
			
			std::vector< int > used_cohorts( cohort_names.size(), 0 ) ;
			std::vector< double > group_sds( group_specs.size(), NA ) ;
			std::vector< std::string > group_names( group_specs.size(), "" ) ;
			std::vector< std::size_t > cohort_group_map( cohort_names.size(), 0 ) ;
			// Fill in the diagonal blocks, one for each group.
			for( std::size_t block_i = 0; block_i < group_specs.size() ; ++block_i ) {
				std::vector< std::string > const elts = split_and_strip( group_specs[block_i], ":" ) ;
				if( elts.size() != 2 ) {
					throw genfile::BadArgumentError(
						"BingwaProcessor::get_group_priors()",
						"model_spec=\"" + group_specs[block_i] + "\"",
						"Group prior spec \"" + group_specs[block_i] + "\" is malformed."
					) ;
				}
				
				// Get the group name...
				std::string const& group_name = elts[0] ;
				group_names[ block_i ] = group_name ;
				// Get the per-group correlation and sd.
				std::vector< boost::tuple< std::vector< std::string >, std::vector< std::string >, std::vector< std::string > > > const covariance_specs = parse_covariance_spec( elts[1] ) ;
				assert( covariance_specs.size() == 1 ) ; // TODO: handle a set of rhos or sds.
				for( std::size_t j = 0; j < covariance_specs.size(); ++j ) {
					boost::tuple< std::vector< std::string >, std::vector< std::string >, std::vector< std::string > > const& covariance_spec = covariance_specs[j] ;
					if( int( covariance_spec.get<eBetweenPopulationCorrelation>().size() ) != m_processor->get_number_of_effect_parameters() ) {
						throw genfile::BadArgumentError(
							"BingwaProcessor::get_group_priors()",
							"model_spec=\"" + model_specs[i] + "\"",
							"In model spec \"" + model_specs[i] + "\", number of correlations and standard deviations specified ("
							+ to_string( between_cohort_rho_string.size() )
							+ ") does not match number of effect size parameters (" + to_string( m_processor->get_number_of_effect_parameters() ) + ")"
						) ;
					}

					double const rho = to_repr< double >( covariance_spec.get<eBetweenPopulationCorrelation>()[0] ) ;
					double const sd = to_repr< double >( covariance_spec.get<eSD>()[0] ) ;
					group_sds[ block_i ] = sd ;
				
					// Get the group members
					std::vector< std::string > const& group_cohorts = get_group_members( group_name, groups ) ;
					std::vector< int > const& cohort_indices = get_cohort_indices( group_cohorts, cohort_names ) ;

					// 
					for( std::size_t k = 0; k < group_cohorts.size(); ++k ) {
						if( used_cohorts[ cohort_indices[k] ] > 0 ) {
							throw genfile::BadArgumentError(
								"BingwaProcessor::get_group_priors()",
								"model_spec=\"" + group_specs[block_i] + "\"",
								"Cohort \"" + group_cohorts[k] + "\" occurs in more than one group."
							) ;
						}
						++used_cohorts[ cohort_indices[k] ] ;
						cohort_group_map[ cohort_indices[k] ] = block_i ;
					}
				
					for( std::size_t cohort_i = 0; cohort_i < group_cohorts.size(); ++cohort_i ) {
						prior( cohort_indices[ cohort_i ], cohort_indices[ cohort_i ] ) = sd*sd ;
						for( std::size_t cohort_j = cohort_i + 1; cohort_j < group_cohorts.size(); ++cohort_j ) {
							prior( cohort_indices[ cohort_i ], cohort_indices[ cohort_j ] ) = sd*sd*rho ;
						}
					}
				}
			}
			// now fill in between-block elements.
			for( std::size_t group1_i = 0; group1_i < group_specs.size() ; ++group1_i ) {
				std::vector< int > const& group1_member_indices = get_cohort_indices( get_group_members( group_names[ group1_i ], groups ), cohort_names ) ;
				for( std::size_t group2_i = group1_i + 1; group2_i < group_specs.size() ; ++group2_i ) {
					std::vector< int > const& group2_member_indices = get_cohort_indices( get_group_members( group_names[ group2_i ], groups ), cohort_names ) ;
					for( std::size_t group1_cohort_i = 0; group1_cohort_i < group1_member_indices.size(); ++group1_cohort_i ) {
						for( std::size_t group2_cohort_i = 0; group2_cohort_i < group2_member_indices.size(); ++group2_cohort_i ) {
							int const i = std::min( group1_member_indices[ group1_cohort_i ], group2_member_indices[ group2_cohort_i ] ) ;
							int const j = std::max( group1_member_indices[ group1_cohort_i ], group2_member_indices[ group2_cohort_i ] ) ;
							prior( i, j ) = between_cohort_rho * group_sds[ group1_i ] * group_sds[ group2_i ] ;
						}
					}
				}
			}
	
			// Finally fill in the lower diagonal
			for( int row = 1; row < prior.rows(); ++row ) {
				for( int col = 0; col < row; ++col ) {
					prior( row, col ) = prior( col, row ) ;
				}
			}
		
			std::string model_name ;
			if( model_names ) {
				model_name = model_names.get()[i] ;
			} else {
				model_name = model_specs[i] ;
			}

			(*result)[ model_name ] = prior ;
		}
	}

	void get_group_specific_priors(
		int const number_of_cohorts,
		GroupDefinition const& groups,
		std::vector< std::string > const& model_specs,
		boost::optional< std::vector< std::string > > model_names,
		std::vector< std::string > const& cohort_names,
		std::map< std::string, Eigen::MatrixXd >* result
	) {
		using namespace genfile::string_utils ;

		if( model_names && model_specs.size() != model_names.get().size() ) {
			throw genfile::BadArgumentError(
				"BingwaProcessor::get_group_specific_priors()",
				"model_names=\"" + join( model_names.get(), " " ) + "\"",
				"Wrong number of group-specific prior names specified (" + to_string( model_names.get().size() ) + ", should be " + to_string( model_specs.size() ) + ")"
			) ;
		}
		
		for( std::size_t i = 0; i < model_specs.size(); ++i ) {
			std::vector< std::string > bits = split_and_strip( model_specs[i], ":" ) ;
			if( bits.size() != 2 ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_group_specific_priors()",
					"model_spec=\"" + model_specs[i] + "\"",
					"Model spec (\"" + model_specs[i] + "\" is malformed."
				) ;
			}
			
			GroupDefinition::const_iterator where = groups.find( bits[0] ) ;
			if( where == groups.end() ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_group_specific_priors()",
					"model_spec=\"" + model_specs[i] + "\"",
					"Group \"" + bits[0] + "\" is not recognised, please define it using the -group option."
				) ;
			}
			std::vector< std::string > const& cohorts = where->second ;

			std::vector< boost::tuple< std::vector< std::string >, std::vector< std::string >, std::vector< std::string > > > const covariance_specs = parse_covariance_spec( bits[1] ) ;
			for( std::size_t j = 0; j < covariance_specs.size(); ++j ) {
				boost::tuple< std::vector< std::string >, std::vector< std::string >, std::vector< std::string > > const& covariance_spec = covariance_specs[j] ;
				if( int( covariance_spec.get<eBetweenPopulationCorrelation>().size() ) != m_processor->get_number_of_effect_parameters() ) {
					throw genfile::BadArgumentError(
						"BingwaProcessor::get_group_specific_priors()",
						"model_spec=\"" + model_specs[i] + "\"",
						"In model spec \"" + model_specs[i] + "\", number of correlations and standard deviations specified ("
						+ to_string( covariance_spec.get<eBetweenPopulationCorrelation>().size() )
						+ ") does not match number of effect size parameters (" + to_string( m_processor->get_number_of_effect_parameters() ) + ")"
					) ;
				}

				double const rho = to_repr< double >( covariance_spec.get<eBetweenPopulationCorrelation>()[0] ) ;
				double const sd = to_repr< double >( covariance_spec.get<eSD>()[0] ) ;
			
				Eigen::MatrixXd prior( number_of_cohorts, number_of_cohorts ) ;
				prior.setZero() ;

				std::vector< int > const& cohort_indices = get_cohort_indices( cohorts, cohort_names ) ;
			
				for( std::size_t j1 = 0; j1 < cohorts.size(); ++j1 ) {
					int const cohort1_j = cohort_indices[j1] ;
					prior( cohort1_j, cohort1_j ) = 1 ;
					for( std::size_t j2 = (j1+1); j2 < cohorts.size(); ++j2 ) {
						int const cohort2_j = cohort_indices[j2] ;
						prior( cohort1_j, cohort2_j ) = prior( cohort2_j, cohort1_j ) = rho ;
					}
				}
			
				prior *= sd * sd ;
			
				std::string model_name ;
				if( model_names ) {
					model_name = model_names.get()[i] ;
				} else {
					model_name = bits[0] + "-specific" ;
				}
			
				(*result)[ model_name + "/rho=" + join( covariance_spec.get<eBetweenPopulationCorrelation>(), " " ) + "/sd=" + join( covariance_spec.get<eSD>(), " " ) ] = prior ;
			}
		}
	}

	std::vector< std::string > parse_rho( std::string const& spec ) {
		using namespace genfile::string_utils ;
		if( spec.substr( 0, 4 ) != "rho=" ) {
			throw genfile::BadArgumentError(
				"BingwaProcessor::parse_rho()",
				"spec=\"" + spec + "\"",
				"Parameter spec \"" + spec + "\" is malformed, should be of the form rho=[rho]/sd=[sd]."
			) ;
		}
		return split_and_strip_discarding_empty_entries( spec.substr( 4, spec.size() ), " " ) ;
	}

	// Given a CovarianceSpec, expand each use of a name from the sets
	//
	template< int element >
	std::vector< CovarianceSpec > expand(
		CovarianceSpec const& covariance_spec,
		std::map< std::string, std::vector< std::string > > const& sets
	) const {
		std::vector< CovarianceSpec > result ;
		// currently we just expand sds according to the sets in sets.
		std::vector< std::string > const& elements = covariance_spec.get< element >() ;
		std::vector< std::string > set_names ;
		std::vector< std::size_t > working_indices ;
		for( std::map< std::string, std::vector< std::string > >::const_iterator i = sets.begin(); i != sets.end(); ++i ) {
			working_indices.push_back( 0 ) ;
			set_names.push_back( i->first ) ;
#if DEBUG_BINGWA
			std::cerr << "expand(): added set name " << set_names.back() << ".\n" ;
#endif
		}
			
		bool complete = false ;
		while( !complete ) {
			// construct the configuration
			std::vector< std::string > these_elements = elements ;
			for( std::size_t i = 0; i < set_names.size(); ++i ) {
				std::string const set_name = set_names[i] ;
				for( std::size_t j = 0; j < elements.size(); ++j ) {
					these_elements[j] = genfile::string_utils::replace_all(
						these_elements[j],
						"[" + set_name + "]",
						sets.find( set_name )->second.at( working_indices[ i ] )
					) ;
#if DEBUG_BINGWA
					std::cerr << "...got \"" << these_elements[j] << "\".\n" ;
#endif
				}
			}
#if DEBUG_BINGWA
			std::cerr << "expand(): these_elements =" ;
			for( std::size_t k = 0; k < these_elements.size(); ++k ) {
				std::cerr << " " << these_elements[k] ;
			}
			std::cerr << "\n" ;
#endif
			CovarianceSpec new_spec( covariance_spec ) ;
			new_spec.get< element >() = these_elements ;
			if( std::find( result.begin(), result.end(), new_spec ) == result.end() ) {
				result.push_back( new_spec ) ;
			}
			
			// move to next configuration of sds and rhos.
			// The rather ungainly code below steps through substitutable values of sds and then of rhos in right-to-left order.
			std::size_t k = 0 ;
			for( ; k < set_names.size(); ++k ) {
				std::size_t j = set_names.size() - k - 1 ;
				std::string const set_name = set_names[j] ;
				if( ++working_indices[j] == sets.find( set_name )->second.size() ) {
					working_indices[j] = 0 ;
				} else {
					break ;
				}
			}
			if( k == set_names.size() ) {
				complete = true ;
			}
		}

		return result ;
	}
	
	template< int element >
	std::vector< CovarianceSpec > expand(
		std::vector< CovarianceSpec > const& covariance_specs,
		std::map< std::string, std::vector< std::string > > const& sets
	) const {
		std::vector< CovarianceSpec > result ;
		for( std::size_t i = 0; i < covariance_specs.size(); ++i ) {
			std::vector< CovarianceSpec > const& this_result = expand< element >( covariance_specs[i], sets ) ;
			result.insert( result.end(), this_result.begin(), this_result.end() ) ;
		}
		return result ;
	}
	
	std::vector< CovarianceSpec > expand_rhos_and_sds(
		CovarianceSpec const& covariance_spec,
		std::map< std::string, std::vector< std::string > > const& sd_sets,
		std::map< std::string, std::vector< std::string > > const& rho_sets
	) const {
		std::vector< CovarianceSpec > result;
		result = expand< eSD >( covariance_spec, sd_sets ) ;
		result = expand< eBetweenPopulationCorrelation >( result, rho_sets ) ;
		return result ;
	}
	
	std::vector< std::string > expand_sd(
		std::string sd,
		std::map< std::string, std::vector< std::string > > const& sd_sets
	) {
		std::vector< std::string > result ;
		std::map< std::string, std::vector< std::string > >::const_iterator where = sd_sets.end() ;
		if( sd.size() > 2 && sd[0] == '[' && sd[ sd.size() -1 ] == ']' ) {
			where = sd_sets.find( sd.substr( 1, sd.size() - 2 )) ;
		}
		if( where == sd_sets.end() ) {
			result.push_back( sd ) ;
		}
		else {
			std::cerr << ">>> Found sd set.\n" ;
			for( std::size_t i = 0; i < where->second.size(); ++i ) {
				result.push_back( where->second[i] ) ;
			}
		}
		return result ;
	}
	
	std::vector< CovarianceSpec > parse_covariance_spec( std::string const& spec ) const {
		using namespace genfile::string_utils ;
		std::vector< std::string > parameters = split_and_strip( spec, "/" ) ;

		int const rho_index = (( parameters.size() > 0 ) && parameters[0].substr(0,4) == "rho=" ) ? 0 : -1 ;
		int const sd_index = (( parameters.size() > (rho_index+1) ) && parameters[ rho_index+1 ].substr(0,3) == "sd=" ) ? (rho_index+1) : -1 ;
		int const cor_index = (( parameters.size() > (sd_index+1) ) && parameters[ sd_index+1 ].substr(0,4) == "cor=" ) ? (sd_index+1) : -1 ;
		bool have_rho = ( rho_index >= 0 ) ;
		bool have_sd = ( sd_index >= 0 ) ;
		bool have_cor = ( cor_index >= 0 ) ;

#if DEBUG_BINGWA
		std::cerr << ( boost::format( "parse_covariance_spec(): spec = %s, rho_index = %d, sd_index = %d, cor_index = %d\n" ) % spec, rho_index, sd_index, cor_index ) ;
#endif		
			
		// Only sd is required.
		if(
			( parameters.size() < 1 || parameters.size() > 3 )
			|| ( !have_sd )
			|| ( parameters.size() == 2 && !( have_cor || have_rho ))
			|| ( parameters.size() == 3 && !( have_cor && have_rho ) )
		) {
			throw genfile::BadArgumentError(
				"BingwaProcessor::parse_covariance_spec()",
				"spec=\"" + spec + "\"",
				"Parameter spec \"" + spec + "\" is malformed, should be of the form "
				"rho=[rho]/sd=[sd]"
				" or rho=[rho]/sd=[sd1 sd1...]/cor=[cor12 cor13...]"
				" or sd=[sd1 sd2...]/cor=[cor11 cor12 cor13...]"
			) ;
		}

		std::vector< std::string > const rhos = ( have_rho ? split_and_strip( parameters[ rho_index ].substr( 4, parameters[ rho_index ].size() ), "," ) : std::vector< std::string >() ) ;
		std::vector< std::string > const sds = ( have_sd ? split_and_strip( parameters[ sd_index ].substr( 3, parameters[ sd_index ].size() ), "," ) : std::vector< std::string >() ) ;
		std::vector< std::string > const cor = ( have_cor ? split_and_strip_discarding_empty_entries( parameters[ cor_index ].substr( 4, parameters[ cor_index ].size() ), "," ) : std::vector< std::string >() ) ;

		if( sds.size() > 1 && !have_cor ) {
			throw genfile::BadArgumentError(
				"BingwaProcessor::parse_covariance_spec()",
				"spec=\"" + spec + "\"",
				"Parameter spec \"" + spec + "\" is malformed, you must supply cor= if there is more than one sd."
			) ;
		}

		std::size_t const expectedNumberOfCorrelations = ( sds.size() * ( sds.size() + 1 ) ) / 2 ;
		if( rhos.size() > 1 ) {
			throw genfile::BadArgumentError(
				"BingwaProcessor::parse_covariance_spec()",
				"spec=\"" + spec + "\"",
				( boost::format( "Wrong number of rhos (%d, should be at most %d)" ) % rhos.size() % 1 ).str()
			) ;
		}
		if( cor.size() != expectedNumberOfCorrelations ) {
			throw genfile::BadArgumentError(
				"BingwaProcessor::parse_covariance_spec()",
				"spec=\"" + spec + "\"",
				( boost::format( "Wrong number of correlations (%d, should be %d)" ) % cor.size() % expectedNumberOfCorrelations ).str()
			) ;
		}
		
		return expand_rhos_and_sds(
			CovarianceSpec( rhos, sds, cor ),
			m_value_sets.at( "sd" ),
			m_value_sets.at( "rho" )
		) ;
	}


	std::vector< CovarianceSpec > expand_sds(
		CovarianceSpec const& covariance_spec,
		std::map< std::string, std::vector< std::string > > const& sd_sets
	) {
		std::vector< CovarianceSpec > result ;
		// currently we just expand sds according to the sets in sd_sets.
		std::vector< std::string > const& sds = covariance_spec.get<eSD>() ;
		std::vector< std::string > sd_set_names ;
		std::vector< std::size_t > working_sd_indices ;
		for( std::map< std::string, std::vector< std::string > >::const_iterator i = sd_sets.begin(); i != sd_sets.end(); ++i ) {
			working_sd_indices.push_back( 0 ) ;
			sd_set_names.push_back( i->first ) ;
		}
			
		bool complete = false ;
		while( !complete ) {
			// construct the configuration
			std::vector< std::string > these_sds = sds ;
			for( std::size_t i = 0; i < sd_set_names.size(); ++i ) {
				std::string const set_name = sd_set_names[i] ;
				for( std::size_t j = 0; j < sds.size(); ++j ) {
					these_sds[j] = genfile::string_utils::replace_all(
						sds[j],
						"[" + sd_set_names[i] + "]",
						sd_sets.find( set_name )->second.at( working_sd_indices[ i ] )
					) ;
				}
			}
			CovarianceSpec new_spec( CovarianceSpec( covariance_spec.get<eBetweenPopulationCorrelation>(), these_sds, covariance_spec.get<eCorrelation>() ) ) ;
			if( std::find( result.begin(), result.end(), new_spec ) == result.end() ) {
				result.push_back( new_spec ) ;
			}
			
			// move to next configuration of sds and rhos.
			// The rather ungainly code below steps through substitutable values of sds and then of rhos in right-to-left order.
			std::size_t k = 0 ;
			for( ; k < sd_set_names.size(); ++k ) {
				std::size_t j = sd_set_names.size() - k - 1 ;
				std::string const set_name = sd_set_names[j] ;
				if( ++working_sd_indices[j] == sd_sets.find( set_name )->second.size() ) {
					working_sd_indices[j] = 0 ;
				} else {
					break ;
				}
			}
			if( k == sd_set_names.size() ) {
				complete = true ;
			}
		}

		return result ;
	}

	std::vector< boost::tuple< std::vector< std::string >, std::vector< std::string >, std::vector< std::string > > > expand_rhos_and_sds(
		std::vector< std::string > rhos,
		std::vector< std::string > sds,
		std::vector< std::string > cor,
		std::map< std::string, std::vector< std::string > > const& sd_sets,
		std::map< std::string, std::vector< std::string > > const& rho_sets
	) const {
		std::vector< boost::tuple< std::vector< std::string >, std::vector< std::string >, std::vector< std::string > > > result ;
		// currently we just expand sds according to the sets in sd_sets.
		std::vector< std::string > working_sds = sds ;
		std::vector< std::string > working_rhos = rhos ;
		std::vector< std::size_t > working_sd_indices( sds.size(), 0 ) ;
		std::vector< std::size_t > working_rho_indices( rhos.size(), 0 ) ;

		bool complete = false ;
		while( !complete ) {
			// construct this configuration
			for( std::size_t j = 0; j < sds.size(); ++j ) {
				if( sds[j].size() > 0 && sds[j][0] == '[' && sds[j][sds[j].size() - 1] == ']' ) {
					std::map< std::string, std::vector< std::string > >::const_iterator where = sd_sets.find( sds[j].substr( 1, sds[j].size() - 2 ) ) ;
					if( where != sd_sets.end() ) {
						working_sds[ j ] = where->second.at( working_sd_indices[j] ) ;
					}
				}

				if( rhos[j].size() > 0 && rhos[j][0] == '[' && rhos[j][rhos[j].size() - 1] == ']' ) {
					std::map< std::string, std::vector< std::string > >::const_iterator where = rho_sets.find( rhos[j].substr( 1, rhos[j].size() - 2 ) ) ;
					if( where != rho_sets.end() ) {
						working_rhos[ j ] = where->second.at( working_rho_indices[j] ) ;
					}
				}
			}
			result.push_back( boost::make_tuple( working_rhos, working_sds, cor ) ) ;
			
			// move to next configuration of sds and rhos.
			// The rather ungainly code below steps through substitutable values of sds and then of rhos in right-to-left order.
			std::size_t k = 0 ;
			for( ; k < sds.size(); ++k ) {
				std::size_t j = sds.size() - k - 1 ;
				if( sds[j].size() > 1 && sds[j][0] == '[' && sds[j][sds[j].size() - 1] == ']' ) {
					std::map< std::string, std::vector< std::string > >::const_iterator where = sd_sets.find( sds[j].substr( 1, sds[j].size() - 2 ) ) ;
					if( where != sd_sets.end() ) {
						working_sd_indices[ j ] = ( working_sd_indices[ j ] + 1 ) % where->second.size() ;
						if( working_sd_indices[ j ] > 0 ) {
							break ;
						}
					}
				}
			}
			if( k == sds.size() ) {
				// No more sds to substitute, so move on to rhos...
				std::size_t k_rho = 0 ;
				for( ; k_rho < rhos.size(); ++k_rho ) {
					std::size_t j = rhos.size() - k_rho - 1 ;
					if( rhos[j].size() > 1 && rhos[j][0] == '[' && rhos[j][rhos[j].size() - 1] == ']' ) {
						std::map< std::string, std::vector< std::string > >::const_iterator where = rho_sets.find( rhos[j].substr( 1, rhos[j].size() - 2 ) ) ;
						if( where != rho_sets.end() ) {
							working_rho_indices[ j ] = ( working_rho_indices[ j ] + 1 ) % where->second.size() ;
							if( working_rho_indices[ j ] > 0 ) {
								break ;
							}
						}
					}
				}
				if( k_rho == rhos.size() ) {
					complete = true ;
				}
			}
		}

		return result ;
	}
	
	std::map< std::string, Eigen::MatrixXd > get_priors( appcontext::OptionProcessor const& options, std::vector< std::string > const& cohort_names ) {
		std::map< std::string, Eigen::MatrixXd > result ;
		using genfile::string_utils::to_string ;
		using genfile::string_utils::to_repr ;
		using genfile::string_utils::split_and_strip_discarding_empty_entries ;
		int const N = options.get_values< std::string >( "-data" ).size() ;
		
		if( options.check( "-define-sd-set" )) {
			m_value_sets[ "sd" ] = parse_value_list( options.get_values< std::string >( "-define-sd-set" )) ;
		} else {
			m_value_sets[ "sd" ] ; // empty map
		}

		if( options.check( "-define-rho-set" )) {
			m_value_sets[ "rho" ] = parse_value_list( options.get_values< std::string >( "-define-rho-set" )) ;
		} else {
			m_value_sets[ "rho" ] ; // empty map
		}

		if( options.check( "-simple-prior" ) ) {
			boost::optional< std::vector< std::string > > model_names ;
			if( options.check( "-simple-prior-name" )) {
				model_names = options.get_values< std::string >( "-simple-prior-name" ) ;
			}
			get_simple_priors(
				N,
				options.get_values< std::string >( "-simple-prior" ),
				model_names,
				&result
			) ;
		}

		if( options.check( "-complex-prior" ) ) {
			get_complex_priors(
				N,
				options.get_values< std::string >( "-complex-prior" ),
				options.get_values< std::string >( "-complex-prior-name" ),
				&result
			) ;
		}

		if( options.check( "-group-specific-prior" ) || options.check( "-group-prior" ) ) {
			GroupDefinition groups ;
			for( std::size_t i = 0; i < cohort_names.size(); ++i ) {
				groups[ cohort_names[i] ] = std::vector< std::string >( 1, cohort_names[i] ) ;
			}
			if( options.check( "-define-group" )) {
				get_groups(
					options.get_values< std::string >( "-define-group"  ),
					cohort_names,
					&groups
				) ;
			}
			
#if DEBUG_BINGWA
			std::cerr << "BingwaApplication:get_priors(): groups are:\n" ;
			GroupDefinition::const_iterator
				i = groups.begin(),
				end_i = groups.end() ;
			for( ; i != end_i; ++i ) {
				std::cerr << "  " << i->first << ":" ;
				for( std::size_t j = 0; j < i->second.size(); ++j ) {
					std::cerr << i->second.at( j ) << " " ;
				}
				std::cerr << "\n" ;
			}
#endif
		
			if( options.check( "-group-prior" )) {
				boost::optional< std::vector< std::string > > model_names ;

				if( options.check( "-group-prior-name" )) {
					model_names = options.get_values< std::string >( "-group-prior-name" ) ;
				}
				get_group_priors( N, groups, options.get_values< std::string >( "-group-prior" ), model_names, cohort_names, &result ) ;
			}

			if( options.check( "-group-specific-prior" )) {
				boost::optional< std::vector< std::string > > model_names ;
				if( options.check( "-group-specific-prior-name" )) {
					model_names = options.get_values< std::string >( "-group-specific-prior-name" ) ;
				}			
				get_group_specific_priors( N, groups, options.get_values< std::string >( "-group-specific-prior" ), model_names, cohort_names, &result ) ;
			}
		}
		
		return result ;
	}

	ValueListSet parse_value_list( std::vector< std::string > const& spec ) const {

		using namespace genfile::string_utils ;
		std::map< std::string, std::vector< std::string > > result ;
		for( std::size_t i = 0; i < spec.size(); ++i ) {
			std::vector< std::string > elts = split_and_strip_discarding_empty_entries( spec[i], "=", " \t" ) ;
			if( elts.size() != 2 ) {
				throw genfile::BadArgumentError(
					"BingwaApplication::parse_value_list()",
					"spec=\"" + spec[i] + "\"",
					"expected specification of the form \"[set name]=[value],[value],[value],...\"."
				) ;
			}
			std::vector< std::string > sds = split_and_strip_discarding_empty_entries( elts[1], ",", " \t" ) ;
			
			// These must be numbers!  Check we can parse them
			for( std::size_t j = 0; j < sds.size(); ++j ) {
				try {
					to_repr< double >( sds[j] ) ;
				}
				catch( StringConversionError const& e ) {
					throw genfile::BadArgumentError(
						"BingwaApplication::parse_value_list()",
						"spec=\"" + spec[i] + "\"",
						"found non-numerical values in value list."
					) ;
				}
			}
			
			result[ elts[0] ] = sds ;
		}
		
#if DEBUG_BINGWA
			std::map< std::string, std::vector< std::string > >::const_iterator
				i = result.begin(),
				end_i = result.end() ;
			std::cerr << "BingwaApplication::parse_value_list(): parsed value lists:\n" ;
			for( ; i != end_i; ++i ) {
				std::cerr << i->first << ": " ;
				for( std::size_t j = 0; j < i->second.size(); ++j ) {
					std::cerr << "\"" << i->second[j] << "\" " ;
				}
				std::cerr << "\n" ;
			}
#endif
		
		return result ;
	}

	void get_groups(
		std::vector< std::string > const& spec,
		std::vector< std::string > const& cohort_names,
		GroupDefinition* result
	) const {
		using namespace genfile::string_utils ;
		assert( result != 0 ) ;
		
		// add default groups
		
		for( std::size_t i = 0; i < spec.size(); ++i ) {
			std::vector< std::string > const bits = split_and_strip( spec[i], "=" ) ;
			if( bits.size() != 2 ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_groups()",
					"spec=\"" + spec[i] + "\"",
					"Group specification \"" + spec[i] + "\" is malformed, should be in the format [name]=group1,group2,..."
				) ;
			}
			if( result->find( bits[0] ) != result->end() ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_groups()",
					"spec=\"" + spec[i] + "\"",
					"Group \"" + bits[0] + "\" is defined more than once (or is the name of a cohort)."
				) ;
			}
			std::vector< std::string > cohorts = split_and_strip( bits[1], "," ) ;
			if( cohorts.size() == 0 ) {
				throw genfile::BadArgumentError(
					"BingwaProcessor::get_groups()",
					"spec=\"" + spec[i] + "\"",
					"Group specification \"" + spec[i] + "\" contains no cohorts"
				) ;
			}
			for( std::size_t i = 0; i < cohorts.size(); ++i ) {
				if( std::find( cohort_names.begin(), cohort_names.end(), cohorts[i] ) == cohort_names.end() ) {
					throw genfile::BadArgumentError(
						"BingwaProcessor::get_groups()",
						"spec=\"" + spec[i] + "\"",
						"Cohort name \"" + cohorts[i] + "\" is not recognised"
					) ;
				}
			}
			
			(*result)[ bits[0] ] = cohorts ;
		}
	}

	Eigen::MatrixXd get_prior_matrix( int const n, double const rho, double const sd ) const {
		return get_correlation_matrix( n, rho ) * sd * sd ;
	}

	Eigen::MatrixXd get_prior_matrix( int const n, std::string const& matrix_spec, double const sd ) const {
		return parse_correlation_matrix( n, matrix_spec ) * sd * sd ;
	}

	Eigen::MatrixXd get_correlation_matrix( int const n, double rho ) const {
		Eigen::MatrixXd result = Eigen::MatrixXd::Identity( n, n ) ;
		for( int i = 0; i < (n-1); ++i ) {
			for( int j = (i+1); j < n; ++j ) {
				result( i, j ) = result( j, i ) = rho ;
			}
		}
		return result ;
	}

	Eigen::MatrixXd parse_correlation_matrix( int const N, std::string const& matrix_spec ) const {
		using namespace genfile::string_utils ;
		Eigen::MatrixXd result = Eigen::MatrixXd::Zero( N, N ) ;
		std::vector< std::string > values = genfile::string_utils::split_and_strip( matrix_spec, "," ) ;
		if( values.size() != (( N * (N+1) ) / 2 ) ) {
			throw genfile::BadArgumentError(
				"BingwaProcessor::parse_correlation_matrix()",
				"matrix_spec=\"" + matrix_spec + "\"",
				"Number of values (" + to_string( values.size() ) + ") is not consistent with the upper triangle of an "
					+ to_string(N) + "x" + to_string(N) + " matrix (should be " + to_string( ( N * (N+1) ) / 2 ) + ")"
			) ;
		}

		{
			int index = 0 ;
			for( int i = 0; i < N; ++i ) {
				for( int j = 0; j < N; ++j ) {
					result( i, j )
						= ( j >= i ) ?
							genfile::string_utils::to_repr< double >( values[ index++ ] )
							:
							result( j, i ) ;
				}
			}
		}
		return result ;
	}

	void summarise_priors(
		std::map< std::string, Eigen::MatrixXd > const& priors,
		std::vector< std::string > const& cohort_names
	) const {
		get_ui_context().logger() << "\n================================================\n" ;
		get_ui_context().logger() << "I will output the following bayesian models:\n\n" ;

		get_ui_context().logger() << std::setprecision( 3 ) ;
		
		std::size_t max_model_name_width = 0 ;
		{
			std::map< std::string, Eigen::MatrixXd >::const_iterator
				prior_i = priors.begin(),
				end_prior_i = priors.end() ;
			for( ; prior_i != end_prior_i; ++prior_i ) {
				max_model_name_width = std::max( max_model_name_width, prior_i->first.size() + 2 ) ;
			}
		}

		EffectParameterNamePack const& parameter_names = m_processor->get_effect_parameter_names() ;
		int const D = parameter_names.size() ;
		int const N = cohort_names.size() ;

		std::vector< std::size_t > column_widths( (D*N), 6 ) ;
		std::size_t max_prefix_width = 6 ;
		for( std::size_t i = 0; i < (D*N); ++i ) {
			std::string const parameter_name = parameter_names.parameter_name( i / N ) ;
			std::string const cohort_name = cohort_names[i % N] ;
			column_widths[i] = 6ul ;//std::max( cohort_name.size() + parameter_name.size() + 2, 6ul ) ;
			max_prefix_width = std::max( max_prefix_width, cohort_name.size() + parameter_name.size() + 7 ) ;
		}
		max_model_name_width = std::max( max_model_name_width, max_prefix_width ) ;

		std::map< std::string, Eigen::MatrixXd >::const_iterator
			prior_i = priors.begin(),
			end_prior_i = priors.end() ;
		std::size_t count = 0 ;
		for( ; prior_i != end_prior_i; ++prior_i, ++count ) {
			get_ui_context().logger() << ( boost::format( "-- (#%.3d) " ) % count ) << prior_i->first << ":\n" ;
			std::string const padding = "  " + std::string( ( max_model_name_width > max_prefix_width ) ? ( max_model_name_width - max_prefix_width ): 0, ' ' ) + "   " ;
			Eigen::MatrixXd const& prior = prior_i->second ;
			assert( prior.rows() == D*N ) ;
			get_ui_context().logger() << "  " << std::setw( max_model_name_width ) << " " << "  " ;
			for( std::size_t i = 0; i < (D*N); ++i ) {
				get_ui_context().logger() << (( i > 0 ) ? " " : "" ) ;
				get_ui_context().logger() << std::setw( column_widths[i] ) << i ;
			}
			get_ui_context().logger() << "\n" ;

			for( int i = 0; i < prior.rows(); ++i ) {
//				get_ui_context().logger() << padding ;
				std::string const parameter_name = parameter_names.parameter_name( i / N ) ;
				std::string const cohort_name = cohort_names[i % N] ;
				get_ui_context().logger()
					<< "  " << std::setw( max_model_name_width )
					<< std::right
					<< ( "(" + cohort_name + ":" + parameter_name + ") " + genfile::string_utils::to_string( i )) << "  "
					<< std::left ;
				for( int j = 0; j < prior.cols(); ++j ) {
					get_ui_context().logger() << (( j > 0 ) ? " " : "" ) ;
					get_ui_context().logger() << std::setw( column_widths[j] ) ;
					if( j < i ) {
						get_ui_context().logger() << " " ;
					} else {
						get_ui_context().logger() << prior(i,j) ;
					}
				}
				get_ui_context().logger() << "\n" ;
			}
			get_ui_context().logger() << "\n";
		}

		get_ui_context().logger() << "================================================\n" ;
	}

	void post_summarise() const {
		std::string const output_file = options().get< std::string >( "-o" ) ;
		std::string const analysis = options().get< std::string >( "-analysis-name" ) ;
		std::string const SQL = "SELECT * FROM SummaryDataView WHERE analysis == '" + analysis + "'" ;
		get_ui_context().logger() << "\n" << globals::program_name << ": Analysis complete.\n" ;
		get_ui_context().logger() << "\n-------------------------\n" ;
		
	}
	
	typedef
		boost::function< void ( genfile::SNPIdentifyingData2 const&, std::string const&, genfile::VariantEntry const& ) > 
		ResultCallback ;

	void load_data() {
		using genfile::string_utils::to_string ;
		using genfile::string_utils::join ;

		std::vector< std::string > cohort_files = options().get_values< std::string >( "-data" ) ;
		for( std::size_t cohort_i = 0; cohort_i < cohort_files.size(); ++cohort_i ) {
			UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading scan results \"" + cohort_files[cohort_i] + "\"" ) ;

			genfile::SNPIdentifyingDataTestConjunction::UniquePtr test( new genfile::SNPIdentifyingDataTestConjunction() ) ;

			if( options().check_if_option_was_supplied_in_group( "SNP inclusion / exclusion options" ) ) {
				genfile::CommonSNPFilter::UniquePtr snp_filter = get_snp_exclusion_filter( cohort_i ) ;
				if( snp_filter.get() ) {
					test->add_subtest( genfile::SNPIdentifyingDataTest::UniquePtr( snp_filter.release() ) ) ;
				}
			}

			boost::optional< std::string > effect_size_column_regex ;
			if( options().check( "-effect-size-column-regex" ) ) {
				std::vector< std::string > values = options().get_values< std::string >( "g-column-regex" ) ;
				if( values.size() != cohort_files.size() ) {
					throw genfile::BadArgumentError(
						"BingwaApplication::load_data()",
						"-effect-size-column-regex=\"" + to_string( genfile::string_utils::join( values, " " ) ) + "\"",
						"Expected " + to_string( cohort_files.size() ) + " values but found " + to_string( values.size() ) + "."
					) ;
				}
				effect_size_column_regex = values[ cohort_i ] ;
			}

			boost::optional< genfile::Chromosome > chromosome_hint ;
			if( options().check( "-assume-chromosome" ) ) {
				chromosome_hint = genfile::Chromosome( options().get< std::string >( "-assume-chromosome" ) ) ;
			}
			FrequentistGenomeWideAssociationResults::UniquePtr results
				= FrequentistGenomeWideAssociationResults::create(
					genfile::wildcard::find_files_by_chromosome( cohort_files[cohort_i] ),
					effect_size_column_regex,
					options().check( "-extra-columns" ) ? options().get_values< std::string >( "-extra-columns" ) : std::vector< std::string >(),
					genfile::SNPIdentifyingDataTest::UniquePtr( test.release() ),
					chromosome_hint,
					SNPTESTResults::SNPResultCallback(),
					progress_context
			) ;
			
			// summarise right now so as to see memory used.
			get_ui_context().logger() << "Cohort " << (cohort_i+1) << " summary: " << results->get_summary() << ".\n" ;
			
			if( cohort_i == 0 ) {
				m_processor->set_effect_parameter_names( results->get_effect_parameter_names() ) ;
			}
			else if( results->get_effect_parameter_names() != m_processor->get_effect_parameter_names() ) {
				throw genfile::MalformedInputError(
					cohort_files[ cohort_i ],
					"Names of effect parameters in cohort " + to_string( cohort_i+1 )
						+ " ("
						+ results->get_effect_parameter_names().get_summary()
						+ ") do not match that in cohort 1 ("
						+ m_processor->get_effect_parameter_names().get_summary()
						+ ")",
					0
				) ;
			}
			m_processor->add_cohort( "cohort_" + to_string( cohort_i+1 ), results ) ;
		}
	}
	
	genfile::CommonSNPFilter::UniquePtr get_snp_exclusion_filter( std::size_t cohort_i ) const {
		genfile::CommonSNPFilter::UniquePtr snp_filter ;
		using genfile::string_utils::join ;
		using genfile::string_utils::split_and_strip_discarding_empty_entries ;

		if( options().check_if_option_was_supplied_in_group( "SNP inclusion / exclusion options" )) {
			snp_filter.reset( new genfile::CommonSNPFilter ) ;

			if( options().check_if_option_was_supplied( "-excl-snpids" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-excl-snpids" ) ;
				BOOST_FOREACH( std::string const& filename, files ) {
					snp_filter->exclude_snps_in_file(
						filename,
						genfile::CommonSNPFilter::SNPIDs
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-incl-snpids" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-incl-snpids" ) ;
				BOOST_FOREACH( std::string const& filename, files ) {
					snp_filter->include_snps_in_file(
						filename,
						genfile::CommonSNPFilter::SNPIDs
					) ;
				}
			}


			if( options().check_if_option_was_supplied( "-excl-rsids" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-excl-rsids" ) ;
				BOOST_FOREACH( std::string const& filename, files ) {
					snp_filter->exclude_snps_in_file(
						filename,
						genfile::CommonSNPFilter::RSIDs
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-incl-rsids" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-incl-rsids" ) ;
				BOOST_FOREACH( std::string const& filename, files ) {
					snp_filter->include_snps_in_file(
						filename,
						genfile::CommonSNPFilter::RSIDs
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-excl-positions" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-excl-positions" ) ;
				BOOST_FOREACH( std::string const& filename, files ) {
					snp_filter->exclude_snps_in_file(
						filename,
						genfile::CommonSNPFilter::Positions
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-incl-positions" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-incl-positions" ) ;
				BOOST_FOREACH( std::string const& filename, files ) {
					snp_filter->include_snps_in_file(
						filename,
						genfile::CommonSNPFilter::Positions
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-excl-snpids-per-cohort" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-excl-snpids-per-cohort" ) ;
				if( files.size() != options().get_values< std::string >( "-data" ).size() ) {
					throw genfile::BadArgumentError( "BingwaApplication::get_snp_exclusion_filter()", "-excl-snpids-per-cohort=\"" + join( files, " " ) + "\"" ) ;
				}
				std::vector< std::string > this_cohort_files = split_and_strip_discarding_empty_entries( files[ cohort_i ], ",", " \t" ) ;
				BOOST_FOREACH( std::string const& filename, this_cohort_files ) {
					snp_filter->exclude_snps_in_file(
						filename,
						genfile::CommonSNPFilter::SNPIDs
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-incl-snpids-per-cohort" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-incl-snpids-per-cohort" ) ;
				if( files.size() != options().get_values< std::string >( "-data" ).size() ) {
					throw genfile::BadArgumentError( "BingwaApplication::get_snp_exclusion_filter()", "-incl-snpids-per-cohort=\"" + join( files, " " ) + "\"" ) ;
				}
				std::vector< std::string > this_cohort_files = split_and_strip_discarding_empty_entries( files[ cohort_i ], ",", " \t" ) ;
				BOOST_FOREACH( std::string const& filename, this_cohort_files ) {
					snp_filter->include_snps_in_file(
						filename,
						genfile::CommonSNPFilter::SNPIDs
					) ;
				}
			}


			if( options().check_if_option_was_supplied( "-excl-rsids-per-cohort" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-excl-rsids-per-cohort" ) ;
				if( files.size() != options().get_values< std::string >( "-data" ).size() ) {
					throw genfile::BadArgumentError( "BingwaApplication::get_snp_exclusion_filter()", "-excl-rsids-per-cohort=\"" + join( files, " " ) + "\"" ) ;
				}
				std::vector< std::string > this_cohort_files = split_and_strip_discarding_empty_entries( files[ cohort_i ], ",", " \t" ) ;
				BOOST_FOREACH( std::string const& filename, this_cohort_files ) {
					snp_filter->exclude_snps_in_file(
						filename,
						genfile::CommonSNPFilter::RSIDs
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-incl-rsids-per-cohort" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-incl-rsids-per-cohort" ) ;
				if( files.size() != options().get_values< std::string >( "-data" ).size() ) {
					throw genfile::BadArgumentError( "BingwaApplication::get_snp_exclusion_filter()", "-incl-rsids-per-cohort=\"" + join( files, " " ) + "\"" ) ;
				}
				std::vector< std::string > this_cohort_files = split_and_strip_discarding_empty_entries( files[ cohort_i ], ",", " \t" ) ;
				BOOST_FOREACH( std::string const& filename, this_cohort_files ) {
					snp_filter->include_snps_in_file(
						filename,
						genfile::CommonSNPFilter::RSIDs
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-excl-snps-per-cohort" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-excl-snps-per-cohort" ) ;
				if( files.size() != options().get_values< std::string >( "-data" ).size() ) {
					throw genfile::BadArgumentError( "BingwaApplication::get_snp_exclusion_filter()", "-excl-snps-per-cohort=\"" + join( files, " " ) + "\"" ) ;
				}
				std::vector< std::string > this_cohort_files = split_and_strip_discarding_empty_entries( files[ cohort_i ], ",", " \t" ) ;
				BOOST_FOREACH( std::string const& filename, this_cohort_files ) {
					genfile::SNPDataSource::UniquePtr source ;
					source.reset(
						new statfile::SNPDataSourceAdapter(
							statfile::BuiltInTypeStatSource::open(
								genfile::wildcard::find_files_by_chromosome( filename )
							)
						)
					) ;

					snp_filter->exclude_snps(
						source->list_snps(),
						genfile::SNPIdentifyingData::CompareFields( options().get_value< std::string >( "-snp-match-fields" ) )
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-incl-snps-per-cohort" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-incl-snps-per-cohort" ) ;
				if( files.size() != options().get_values< std::string >( "-data" ).size() ) {
					throw genfile::BadArgumentError( "BingwaApplication::get_snp_exclusion_filter()", "-incl-snps-per-cohort=\"" + join( files, " " ) + "\"" ) ;
				}
				std::vector< std::string > this_cohort_files = split_and_strip_discarding_empty_entries( files[ cohort_i ], ",", " \t" ) ;
				BOOST_FOREACH( std::string const& filename, this_cohort_files ) {
					genfile::SNPDataSource::UniquePtr source ;
					source.reset(
						new statfile::SNPDataSourceAdapter(
							statfile::BuiltInTypeStatSource::open(
								genfile::wildcard::find_files_by_chromosome( filename )
							)
						)
					) ;

					snp_filter->include_snps(
						source->list_snps(),
						genfile::SNPIdentifyingData::CompareFields( options().get_value< std::string >( "-snp-match-fields" ) )
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-excl-snps-matching" )) {
				std::string const it = options().get< std::string > ( "-excl-snps-matching" ) ;
				std::vector< std::string > specs = genfile::string_utils::split_and_strip_discarding_empty_entries( it, ",", " \t" ) ;
				BOOST_FOREACH( std::string const& spec, specs ) {
					snp_filter->exclude_snps_matching(
						spec
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-incl-snps-matching" )) {
				std::string const it = options().get< std::string > ( "-incl-snps-matching" ) ;
				std::vector< std::string > specs = genfile::string_utils::split_and_strip_discarding_empty_entries( it, ",", " \t" ) ;
				BOOST_FOREACH( std::string const& spec, specs ) {
					snp_filter->include_snps_matching(
						spec
					) ;
				}
			}
			
			if( options().check_if_option_was_supplied( "-incl-range" )) {
				std::vector< std::string > const specs = options().get_values< std::string >( "-incl-range" ) ;
				for ( std::size_t i = 0; i < specs.size(); ++i ) {
					snp_filter->include_snps_in_range(
						genfile::GenomePositionRange::parse( specs[i] )
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-excl-range" )) {
				std::vector< std::string > specs = options().get_values< std::string >( "-excl-range" ) ;
				for ( std::size_t i = 0; i < specs.size(); ++i ) {
					snp_filter->exclude_snps_in_range(
						genfile::GenomePositionRange::parse( specs[i] )
					) ;
				}
			}
		}
		return snp_filter ;
	}
	
private:
	std::map< std::string, ValueListSet > m_value_sets ;
	BingwaProcessor::UniquePtr m_processor ;
} ;

int main( int argc, char **argv ) {
	try {
		BingwaApplication app( argc, argv ) ;	
		app.run() ;
		app.post_summarise() ;
	}
	catch( appcontext::HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}
