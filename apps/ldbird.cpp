
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <iostream>
#include <limits>
#include <typeinfo>
#include <fstream>
#include <unordered_map>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/timer/timer.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include "config/qctool_version_autogenerated.hpp"

#include "appcontext/CmdLineOptionProcessor.hpp"
#include "appcontext/ApplicationContext.hpp"

#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPDataSourceChain.hpp"
#include "genfile/ThreshholdingSNPDataSource.hpp"
#include "genfile/ToGP.hpp"
#include "genfile/CommonSNPFilter.hpp"
#include "genfile/VariantIdentifyingDataFilteringSNPDataSource.hpp"
#include "genfile/CrossCohortCovariateValueMapping.hpp"
#include "genfile/SampleFilter.hpp"
#include "genfile/SampleFilteringSNPDataSource.hpp"
#include "genfile/SampleFilteringCohortIndividualSource.hpp"
#include "genfile/CompoundSampleFilter.hpp"
#include "genfile/SampleFilterNegation.hpp"
#include "genfile/VariableInSetSampleFilter.hpp"
#include "genfile/GPThresholdingGTSetter.hpp"
#include "genfile/db/Error.hpp"

#include "components/SNPSummaryComponent/InfoComputation.hpp"

#include "metro/constants.hpp"
#include "metro/regression/Design.hpp"
#include "metro/SampleRange.hpp"
#include "metro/intersect_ranges.hpp"
#include "metro/regression/BinomialLogistic.hpp"
#include "metro/regression/ThreadedLogLikelihood.hpp"
#include "metro/regression/IndependentNormalWeightedLogLikelihood.hpp"
#include "metro/regression/IndependentLogFWeightedLogLikelihood.hpp"
#include "metro/regression/LogPosteriorDensity.hpp"
#include "metro/ValueStabilisesStoppingCondition.hpp"
#include "metro/Snptest25StoppingCondition.hpp"
#include "metro/maximisation.hpp"
#include "metro/fit_model.hpp"
#include "metro/CholeskyStepper.hpp"

#include "qcdb/MultiVariantStorage.hpp"
#include "qcdb/FlatTableDBOutputter.hpp"

#include <Eigen/Core>
#include <Eigen/QR>
#include "Eigen/Eigenvalues"

// #define DEBUG 1

namespace globals {
	std::string const program_name = "ldbird" ;
	std::string const program_version = qctool_version ;
	std::string const program_revision =  std::string( qctool_revision ).substr( 0, 7 ) ;
}

namespace {
	std::vector< std::string > collect_unique_ids( std::vector< std::string > const& ids_or_filenames ) {
		std::vector< std::string > result ;
		for( auto elt: ids_or_filenames ) {
			if( boost::filesystem::exists( elt )) {
				std::ifstream f( elt ) ;
				std::copy(
					std::istream_iterator< std::string >( f ),
					std::istream_iterator< std::string >(),
					std::back_inserter< std::vector< std::string > >( result )
				) ;
			} else {
				result.push_back( elt ) ;
			}
		}
		// now sort and uniqueify them...
		std::sort( result.begin(), result.end() ) ;
		std::vector< std::string >::const_iterator newBack = std::unique( result.begin(), result.end() ) ;
		result.resize( newBack - result.begin() ) ;
		return result ;
	}

	int extract_mapped_categorical_value(
		std::size_t sample_index,
		genfile::CohortIndividualSource const& samples,
		std::size_t column_index,
		genfile::CrossCohortCovariateValueMapping const& mapping
	) {
		int value = -1 ;
		genfile::VariantEntry entry = samples.get_entry( sample_index, column_index ) ;
		if( !entry.is_missing() ) {
			value = mapping.get_mapped_value( entry ).as< int >() ;
		}
		return value ;
	}

	std::string get_unmapped_level(
		genfile::CrossCohortCovariateValueMapping const& mapping,
		int level
	) {
		genfile::VariantEntry entry = mapping.get_unmapped_value( level ) ;
		return entry.as< std::string >() ;
	}
}

struct Visitor: public boost::noncopyable {
public:
	typedef boost::function<
		void (
			std::vector< int > const& changed,
			std::vector< genfile::VariantIdentifyingData > const& variants,
			std::vector< genfile::VariantDataReader::SharedPtr > const& readers
	) > Callback ;
public:
	~Visitor() {} ;
	virtual bool step( Callback callback ) = 0 ;
	virtual boost::optional< std::size_t > count() const = 0 ;
} ;

struct CartesianProduct: public Visitor {
public:
	CartesianProduct( bool lower_triangle = false ):
		m_lower_triangle( lower_triangle )
	{}

	void add_source(
		std::string const& name,
		genfile::SNPDataSource* source
	) {
		m_names.push_back( name ) ;
		m_sources.push_back( source ) ;
	}
	
	bool step( Callback callback ) {
		if( m_sources.empty() ) {
			return false ;
		}
		if( m_changed.empty() ) {
			// First call.  Populate our data structures.
			m_changed.assign( m_sources.size(), 1 ) ;
			m_variants.resize( m_sources.size() ) ;
			assert( m_readers.size() == 0 ) ;
			for( std::size_t i = 0; i < m_sources.size(); ++i ) {
				if( !m_sources[i]->get_snp_identifying_data( &m_variants[i] )) {
					return false ;
				} ;
				m_readers.push_back( m_sources[i]->read_variant_data() ) ;
			}
		} else {
			// Move to the next variant.
			// We step to the next variant in each source, starting with the last.
			// If a source is exhausted we reset it to the start and recurse to
			// the preceding source.
			std::size_t i = m_sources.size() - 1 ;
			for( ; true; --i ) {
				if( !m_sources[i]->get_snp_identifying_data( &m_variants[i] )) {
					if( i == 0 ) {
						return false ;
					} else {
						reset_source(i) ;
						m_sources[i]->get_snp_identifying_data( &m_variants[i] ) ;
					}
				} else {
					break ;
				}
			}
			m_changed.assign( m_sources.size(), 0 ) ;
			for( std::size_t j = i; j < m_sources.size(); ++j ) {
				m_changed[j] = 1 ;
				m_readers[j] = m_sources[j]->read_variant_data() ;
			}
		}
		callback(
			m_changed,
			m_variants,
			m_readers
		) ;
		return true ;
	}
	
	boost::optional< std::size_t > count() const {
		boost::optional< std::size_t > result ;
		if( m_sources.size() > 0 && (result = m_sources[0]->size() )) {
			if( m_lower_triangle ) {
				assert( m_sources.size() == 2 ) ;
				result = ((*result) * (*result + 1) / 2) ;
			} else {
				for( std::size_t i = 1; i < m_sources.size(); ++i ) {
					boost::optional< std::size_t > source_size = m_sources[i]->size() ;
					if( !source_size ) {
						return boost::optional< std::size_t >() ;
					} else {
						(*result) *= (*source_size) ;
					}
				}
			}
		}
		return result ;
	}
	
private:
	std::vector< std::string > m_names ;
	bool const m_lower_triangle ;
	std::vector< genfile::SNPDataSource* > m_sources ;
	std::vector< genfile::VariantIdentifyingData > m_variants ;
	std::vector< genfile::VariantDataReader::SharedPtr > m_readers ;
	std::vector< int > m_changed ;
	
private:
	
	void reset_source( std::size_t i ) {
		m_sources[i]->reset_to_start() ;
		if( m_lower_triangle ) {
			assert( m_sources.size() == 2 ) ;
			if( i == 1 ) {
				genfile::VariantIdentifyingData variant ;
				while( m_sources[1]->number_of_snps_read() < m_sources[0]->number_of_snps_read() ) {
					m_sources[1]->get_snp_identifying_data( &variant ) ;
					m_sources[1]->ignore_snp_probability_data() ;
				}
			}
		}
	}
} ;

struct CorrelatorOptions: public appcontext::CmdLineOptionProcessor
{
public:
	std::string get_program_name() const { return globals::program_name ; }

	void declare_options( appcontext::OptionProcessor& options ) {
		// Meta-options
		options.set_help_option( "-help" ) ;
		options.set_spec_option( "-spec" ) ;

		// File options
		options.declare_group( "Input file options" ) ;
		options[ "-g1" ]
			.set_description( 	"Path to first genotype file."
								"The given filename may contain the wildcard character '#', which expands to match a"
								"one- or two-character chromosome identifier." )
			.set_takes_values( 1 )
			.set_minimum_multiplicity( 1 )
			.set_maximum_multiplicity( 1 ) ;
		options[ "-g1-incl-rsids" ]
			.set_description( "Exclude all SNPs whose RSID is not in the given file(s) from the analysis.")
			.set_takes_values_until_next_option()
			.set_maximum_multiplicity( 100 ) ;
		options[ "-g1-incl-range" ]
			.set_description( "Specify a range of SNPs (or comma-separated list of ranges of SNPs) to operate on. "
				"Each range should be in the format CC:xxxx-yyyy where CC is the chromosome and xxxx and yyyy are the "
				"start and end coordinates, or just xxxx-yyyy which matches that range from all chromosomes. "
				"You can also omit either of xxxx or yyyy to get all SNPs from the start or to the end of a chromosome." )
			.set_takes_values_until_next_option() ;

		options.option_implies_option( "-g1-incl-range", "-g1" ) ;
		options.option_implies_option( "-g1-incl-rsids", "-g1" ) ;

		options[ "-g2" ]
			.set_description( 	"Path to second genotype file.  If not given the first genotype file will be used."
								"The given filename may contain the wildcard character '#', which expands to match a"
								"one- or two-character chromosome identifier." )
			.set_takes_values( 1 )
			.set_minimum_multiplicity( 0 )
			.set_maximum_multiplicity( 1 ) ;

		options[ "-g2-incl-rsids" ]
			.set_description( "Exclude all SNPs whose RSID is not in the given file(s) from the analysis.")
			.set_takes_values_until_next_option()
			.set_maximum_multiplicity( 100 ) ;
		options[ "-g2-incl-range" ]
			.set_description( "Specify a range of SNPs (or comma-separated list of ranges of SNPs) to operate on. "
				"Each range should be in the format CC:xxxx-yyyy where CC is the chromosome and xxxx and yyyy are the "
				"start and end coordinates, or just xxxx-yyyy which matches that range from all chromosomes. "
				"You can also omit either of xxxx or yyyy to get all SNPs from the start or to the end of a chromosome." )
			.set_takes_values_until_next_option() ;

		options.option_implies_option( "-g2-incl-range", "-g2" ) ;
		options.option_implies_option( "-g2-incl-rsids", "-g2" ) ;

		options[ "-s" ]
			.set_description( "Path of sample file" )
			.set_takes_values( 1 )
			.set_minimum_multiplicity( 1 )
			.set_maximum_multiplicity( 1 ) ;

		options.declare_group( "Output file options" ) ;
		options[ "-o" ]
			.set_description( "Output file" )
			.set_takes_values( 1 )
			.set_minimum_multiplicity( 1 )
			.set_maximum_multiplicity( 1 ) ;
		
		options.declare_group( "Model options" ) ;
		options[ "-incl-samples"]
			.set_description( "Filter out samples whose sample ID does not lie in the given file(s).")
			.set_takes_values( 1 )
			.set_minimum_multiplicity( 0 )
			.set_maximum_multiplicity( 100 ) ;
		options[ "-excl-samples"]
			.set_description( "Filter out samples whose sample ID lies in the given file.")
			.set_takes_values( 1 )
			.set_minimum_multiplicity( 0 )
			.set_maximum_multiplicity( 100 ) ;
		options[ "-incl-samples-where"]
			.set_description( "Include samples by specifying conditions on the values of columns in the sample file.")
			.set_takes_values( 1 )
			.set_minimum_multiplicity( 0 )
			.set_maximum_multiplicity( 100 ) ;
		options[ "-excl-samples-where"]
			.set_description( "Exclude samples by specifying conditions on the values of columns in the sample file.")
			.set_takes_single_value()
			.set_minimum_multiplicity( 0 )
			.set_maximum_multiplicity( 100 ) ;
		options[ "-minimum-r2" ]
			.set_description( "Do not output results where r^2 is lower than this threshold." )
			.set_takes_single_value()
			.set_default_value( 0.05 ) ;
		options[ "-assume-haploid" ]
			.set_description( "Convert all data to haploid calls.  This converts homozygous calls"
				" to haploid calls, and treats any heterozygous calls as missing.")
		;
		options[ "-prior-weight" ]
			.set_description( "Specify a prior weight for samples.  "
				" This augments data with one of each of the four possible haploid genotype combinations "
				" for each pairwise comparison, with the total weight given."
				" The default value (1) indicates each haplotype is given 1/4 weight."
			)
			.set_takes_single_value()
			.set_default_value( "1.0" )
		;
		options.declare_group( "Miscellaneous options" ) ;
		options[ "-analysis-name" ]
			.set_description( "Specify a name to label results from this analysis with." )
			.set_takes_single_value()
			.set_default_value( "qctool analysis" ) ;
		options[ "-analysis-chunk" ]
			.set_description( "Specify a name denoting the current genomic region or chunk on which this is run.  This is intended for use in parallel environments." )
			.set_takes_single_value()
			.set_default_value( genfile::MissingValue() ) ;
		options[ "-table-prefix" ]
			.set_description( "Specify a prefix to add to tables.  They will be called <prefix>Frequency and <prefix>R." )
			.set_takes_single_value()
			.set_default_value( "" ) ;
		options[ "-threshold" ]
			.set_description( "The threshold to apply to genotype probabilities, if necessary, to make calls." )
			.set_takes_single_value()
			.set_default_value( 0.9 ) ;
		options.declare_group( "Miscellaneous options" ) ;
		options[ "-debug" ]
			.set_description( "Output debugging information." ) ;
	}
} ;


namespace {
	struct ProbSetter: public genfile::VariantDataReader::PerSampleSetter {
		typedef std::vector< metro::SampleRange > SampleRanges ;

		ProbSetter(
			Eigen::MatrixXd* genotypes,
			Eigen::VectorXi* ploidy,
			SampleRanges* nonmissing_samples
		):
			m_genotypes( genotypes ),
			m_ploidy( ploidy ),
			m_nonmissing_samples( nonmissing_samples ),
			m_sample_i(0)
		{
			assert( genotypes != 0 && nonmissing_samples != 0 && ploidy != 0 ) ;
		}

		void initialise( std::size_t nSamples, std::size_t nAlleles ) {
			if( nAlleles != 2 ) {
				throw genfile::BadArgumentError(
					"ProbSetter::initialise()",
					"nAlleles=" + genfile::string_utils::to_string( nAlleles ),
					"I only support biallelic variants"
				) ;
			}
			m_genotypes->resize( nSamples, 3 ) ;
			m_genotypes->setZero() ;
			m_ploidy->resize( nSamples ) ;
			m_ploidy->setZero() ;
			
			m_nonmissing_samples->clear() ;
			m_last_nonmissing_sample_i = 0 ;
		}

		bool set_sample( std::size_t n ) {
			m_sample_i = n ;
			return true ;
		}

		void set_number_of_entries(
			uint32_t ploidy, std::size_t n,
			genfile::OrderType const order_type,
			genfile::ValueType const value_type
		) {
			if( value_type != genfile::eProbability ) {
				throw genfile::BadArgumentError(
					"Callsetter::set_number_of_entries()",
					"value_type=" + genfile::string_utils::to_string( value_type ),
					"Expected genotype call probabilities (e.g. GP field)."
				) ;
			}
			if( order_type != genfile::ePerUnorderedGenotype ) {
				throw genfile::BadArgumentError(
					"Callsetter::set_number_of_entries()",
					"order_type=" + genfile::string_utils::to_string( order_type ),
					"Expected genotype call probabilities (e.g. GP field)."
				) ;
			}
			if( ploidy != 2 ) {
				throw genfile::BadArgumentError(
					"Callsetter::set_number_of_entries()",
					"order_type=" + genfile::string_utils::to_string( order_type ),
					"Expected a diploid sample."
				) ;
			}
			(*m_ploidy)(m_sample_i) = ploidy ;
		}

		void set_value( std::size_t entry_i, genfile::MissingValue const value ) {
			// This sample has missing data, end our current sample range here.
			if( m_sample_i > m_last_nonmissing_sample_i ) {
				m_nonmissing_samples->push_back(
					metro::SampleRange( m_last_nonmissing_sample_i, m_sample_i )
				) ;
			}
			// Skip this sample for next range.
			m_last_nonmissing_sample_i = m_sample_i + 1 ;
		}

		void set_value( std::size_t entry_i, Integer const value ) {
			assert(0) ; // expected a probabilty
		}

		void set_value( std::size_t entry_i, double const value ) {
			if( m_sample_i >= m_last_nonmissing_sample_i ) {
				(*m_genotypes)(m_sample_i, entry_i) = value ;
			}
		}

		void finalise() {
			// Add the final sample range if needed.
			if( (m_sample_i+1) > m_last_nonmissing_sample_i ) {
				m_nonmissing_samples->push_back(
					metro::SampleRange( m_last_nonmissing_sample_i, m_sample_i+1 )
				) ;
			}
			
			// For probability data, we need to handle the case of missingness encoded as zero probabilities.
			// Do that now by inspecting the genotypes.
			{
				std::vector< metro::SampleRange > ranges ;
				std::size_t last_nonmissing_sample_i = 0 ;
				std::size_t i = 0 ;
				for( ; i < m_genotypes->rows(); ++i ) {
					if( m_genotypes->row(i).sum() == 0.0 ) {
						if( i > last_nonmissing_sample_i ) {
							ranges.push_back( metro::SampleRange( last_nonmissing_sample_i, i )) ;
						}
						last_nonmissing_sample_i = i+1 ;
					}
				}
				if( i > last_nonmissing_sample_i ) {
					ranges.push_back( metro::SampleRange( last_nonmissing_sample_i, i )) ;
				}
				
				*m_nonmissing_samples = metro::impl::intersect_ranges( *m_nonmissing_samples, ranges ) ;
			}
		}
		
	private:
		Eigen::MatrixXd* m_genotypes ;
		Eigen::VectorXi* m_ploidy ;
		SampleRanges* m_nonmissing_samples ;
		std::size_t m_last_nonmissing_sample_i ;
		std::size_t m_sample_i ;
	} ;
}

struct CorrelatorApplication: public appcontext::ApplicationContext
{
public:
	
	typedef Eigen::MatrixXd Matrix ;
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::VectorXi IntegerVector ;
	typedef Eigen::RowVectorXd RowVector ;
	
public:
	CorrelatorApplication( int argc, char** argv ):
		appcontext::ApplicationContext(
			globals::program_name,
			globals::program_version + ", revision " + globals::program_revision,
			std::auto_ptr< appcontext::OptionProcessor >( new CorrelatorOptions ),
			argc,
			argv,
			"-log"
		)
	{
		process() ;
	}
	
private:
	void process() {
		try {
			unsafe_process() ;
		}
		catch( genfile::InputError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): " << e.format_message() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( genfile::FileNotFoundError const& e ) {
			get_ui_context().logger() << "\nError: No file matching \"" << e.filespec() << "\" could be found.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( genfile::db::Error const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << ") with the following statement: \""
				<< e.sql()
				<< "\".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
	}
	
	void unsafe_process() {
		using genfile::string_utils::to_string ;
		
		genfile::SampleFilterConjunction::UniquePtr sample_filter = get_sample_filter() ;

		genfile::CohortIndividualSource::UniquePtr
			samples = genfile::CohortIndividualSource::create( options().get< std::string >( "-s" ) ) ;

		std::set< std::size_t > const excluded_samples = compute_excluded_samples( *sample_filter, *samples ) ;
		if( excluded_samples.size() > 0 ) {
			samples.reset(
				genfile::SampleFilteringCohortIndividualSource::create( samples, excluded_samples ).release()
			) ;
		}

		genfile::SNPDataSource::UniquePtr g1 = open_genotype_data_sources(
			options().get< std::string >( "-g1" ),
			get_variant_filter(
				"-g1-incl-range",
				"-g1-incl-rsids"
			),
			excluded_samples
		) ;

		genfile::SNPDataSource::UniquePtr g2 ;
		if( options().check( "-g2" )) {
			// Two files given, we'll do a full carThis mechanism chooses a full cartesian product or a lower triangle implementation.
			// TODO: make this more obvious / cleaner.
			g2 = open_genotype_data_sources(
				options().get< std::string >( "-g2" ),
				get_variant_filter(
					"-g2-incl-range",
					"-g2-incl-rsids"
				),
				excluded_samples
			) ;
		} else {
			g2 = open_genotype_data_sources(
				options().get< std::string >( "-g1" ),
				get_variant_filter(
					"-g1-incl-range",
					"-g1-incl-rsids"
				),
				excluded_samples
			) ;
		}
		
		std::size_t const N = samples->size() ;
		if( g1->number_of_samples() != N ) {
			throw genfile::BadArgumentError(
				"CorrelatorApplication::unsafe_process()",
				"-g1 \"" + options().get< std::string >( "-g1" ) + "\"",
				"Wrong number of samples (" + to_string(g1->number_of_samples()) + ", expected " + to_string(N) + ")"
			) ;
		}
		if( g2->number_of_samples() != N ) {
			throw genfile::BadArgumentError(
				"CorrelatorApplication::unsafe_process()",
				"-g2 \"" + options().get< std::string >( "-g2" ) + "\"",
				"Wrong number of samples (" + to_string(g2->number_of_samples()) + ", expected " + to_string(N) + ")"
			) ;
		}
		
		write_preamble( *g1, *g2, *samples ) ;
		
		{
			std::string const fileSpec = options().get< std::string > ( "-o" ) ;
			std::string const tablePrefix = options().get< std::string > ( "-table-prefix" ) ;
			qcdb::Storage::UniquePtr frequencyStorage = qcdb::Storage::create(
				fileSpec + ":" + tablePrefix + "Frequency",
				options().get< std::string > ( "-analysis-name" ),
				options().get< std::string > ( "-analysis-chunk" ),
				get_application_metadata()
			) ;
			//frequencyStorage->set_variant_names( std::vector< std::string >({ "g1", "g2" })) ;

			qcdb::MultiVariantStorage::UniquePtr correlationStorage = qcdb::MultiVariantStorage::create(
				fileSpec + ":" + tablePrefix + "R+without rowid",
				2,
				options().get< std::string > ( "-analysis-name" ),
				options().get< std::string > ( "-analysis-chunk" ),
				get_application_metadata()
			) ;
			correlationStorage->set_variant_names( std::vector< std::string >({ "g1", "g2" })) ;

			// This mechanism chooses a full cartesian product or a lower triangle implementation
			// TODO: make this more obvious / cleaner.
			CartesianProduct visitor( !options().check( "-g2" ) ) ;
			visitor.add_source( "g1", g1.get() ) ;
			visitor.add_source( "g2", g2.get() ) ;

			run( visitor, *samples, *frequencyStorage, *correlationStorage ) ;

			frequencyStorage->finalise() ;
			correlationStorage->finalise() ;
		}
	}
	
	genfile::CommonSNPFilter::UniquePtr get_variant_filter(
		std::string const& incl_range_option,
		std::string const& incl_ids_option
	) const {
		genfile::CommonSNPFilter::UniquePtr result ;
		if( options().check( incl_range_option ) || options().check( incl_ids_option ) ) {
			result.reset( new genfile::CommonSNPFilter() ) ;
			if( options().check( incl_range_option )) {
				std::vector< std::string > specs = collect_unique_ids( options().get_values< std::string >( incl_range_option ) ) ;
				for ( std::size_t i = 0; i < specs.size(); ++i ) {
					result->include_snps_in_range(
						genfile::GenomePositionRange::parse( specs[i] )
					) ;
				}
			}
			if( options().check( incl_ids_option )) {
				std::vector< std::string > files = collect_unique_ids( options().get_values< std::string > ( incl_ids_option ) ) ;
				result->include_snps_in_set(
					std::set< std::string >( files.begin(), files.end() ),
					genfile::CommonSNPFilter::RSIDs
				) ;
			}
		}
		return result ;
	}
	
	genfile::SampleFilterConjunction::UniquePtr get_sample_filter() const {
		genfile::SampleFilterConjunction::UniquePtr filter( new genfile::SampleFilterConjunction() ) ;
		genfile::SampleFilterDisjunction::UniquePtr sample_inclusion_filter( new genfile::SampleFilterDisjunction() ) ;

		if( options().check_if_option_was_supplied( "-incl-samples" ) ) {
			genfile::SampleFilter::UniquePtr id1_filter = get_sample_id_filter(
				options().get_values< std::string >( "-incl-samples" )
			) ;
			sample_inclusion_filter->add_clause( genfile::SampleFilter::UniquePtr( id1_filter.release() ) ) ;
		}

		if( options().check_if_option_was_supplied( "-excl-samples" ) ) {
			genfile::SampleFilter::UniquePtr id1_filter = get_sample_id_filter(
				options().get_values< std::string >( "-excl-samples" )
			) ;
			filter->add_clause(
				genfile::SampleFilter::UniquePtr(
					new genfile::SampleFilterNegation(
						id1_filter
					)
				)
			) ;
		}

		if( options().check_if_option_was_supplied( "-incl-samples-where" ) ) {
			std::vector< std::string > conditions = options().get_values< std::string >( "-incl-samples-where" ) ;
			// we OR together different WHERE clauses.
			for( std::size_t i = 0; i < conditions.size(); ++i ) {
				sample_inclusion_filter->add_clause( genfile::SampleFilter::create( conditions[i] )) ;
			}
		}

		if( options().check_if_option_was_supplied( "-excl-samples-where" ) ) {
			std::vector< std::string > conditions = options().get_values< std::string >( "-excl-samples-where" ) ;
			// we AND together different WHERE clauses.
			genfile::SampleFilterDisjunction::UniquePtr where( new genfile::SampleFilterDisjunction() ) ;
			for( std::size_t i = 0; i < conditions.size(); ++i ) {
				where->add_clause( genfile::SampleFilter::create( conditions[i] )) ;
			}
			filter->add_clause(
				genfile::SampleFilter::UniquePtr(
					new genfile::SampleFilterNegation(
						genfile::SampleFilter::UniquePtr( where.release() )
					)
				)
			) ;
		}
		if( sample_inclusion_filter->number_of_clauses() > 0 ) {
			filter->add_clause( genfile::SampleFilter::UniquePtr( sample_inclusion_filter.release() ) ) ;
		}

		return filter ;
	}

	std::set< std::size_t > compute_excluded_samples( genfile::SampleFilter const& filter, genfile::CohortIndividualSource const& samples ) const {
		std::vector< std::size_t > included_samples ;
		void(std::vector<std::size_t>::*push_back)(std::size_t const&) = &std::vector<std::size_t>::push_back ;
		filter.test(
			samples,
			boost::bind(
				push_back,
				&included_samples,
				_1
			)
		) ;
		std::set< std::size_t > excluded_samples(
			boost::counting_iterator< std::size_t >(0),
			boost::counting_iterator< std::size_t >( samples.get_number_of_individuals() )
		) ;
		for( std::size_t i = 0; i < included_samples.size(); ++i ) {
			excluded_samples.erase( included_samples[i] ) ;
		}
		return excluded_samples ;
	}

	genfile::SampleFilter::UniquePtr get_sample_id_filter( std::vector< std::string > const& filenames ) const {
		genfile::VariableInSetSampleFilter::UniquePtr result( new genfile::VariableInSetSampleFilter( "ID_1" )) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			std::string elt ;
			std::ifstream is( filenames[i].c_str() ) ;
			while( is >> elt ) {
				result->add_level( elt ) ;
			}
		}
		return genfile::SampleFilter::UniquePtr( result.release() ) ;
	}
	
	genfile::SNPDataSource::UniquePtr open_genotype_data_sources(
		std::string const& filename,
		genfile::CommonSNPFilter::UniquePtr filter,
		std::set< std::size_t > const& excluded_samples,
		bool threshold = false
	) {
		std::vector< genfile::wildcard::FilenameMatch > filenames
			= genfile::wildcard::find_files_by_chromosome(
				filename,
				genfile::wildcard::eALL_CHROMOSOMES
			) ;

		genfile::SNPDataSource::UniquePtr source( 
			genfile::SNPDataSourceChain::create(
				filenames
			).release()
		) ;
		if( filter.get() ) {
			source.reset(
				genfile::VariantIdentifyingDataFilteringSNPDataSource::create(
					source,
					genfile::VariantIdentifyingDataTest::UniquePtr( filter.release() )
				).release()
			) ;
		}
		if( threshold ) {
			source.reset(
				new genfile::ThreshholdingSNPDataSource( source, 0.9 )
			) ;
		}
		if( excluded_samples.size() > 0 ) {
			source.reset(
				new genfile::SampleFilteringSNPDataSource( source, excluded_samples )
			) ;
		}
		return source ;
	}

	void write_preamble(
		genfile::SNPDataSource const& host,
		genfile::SNPDataSource const& para,
		genfile::CohortIndividualSource const& samples
	) {
		get_ui_context().logger() << "Loaded data for " << samples.size() << " samples.\n" ;
		get_ui_context().logger() << "  First dataset:\n" << host.get_summary() << "\n" ;
		get_ui_context().logger() << " Second dataset:\n" << para.get_summary() << "\n" ;
	}


	void run(
		Visitor& visitor,
		genfile::CohortIndividualSource const& samples,
		qcdb::Storage& frequencyOutput,
		qcdb::MultiVariantStorage& correlationOutput
	) {
		appcontext::UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Testing" ) ;
		for(
			std::size_t count = 0;
			visitor.step(
				[this,&frequencyOutput,&correlationOutput](
					std::vector< int > const& changed,
					std::vector< genfile::VariantIdentifyingData > const& variants,
					std::vector< genfile::VariantDataReader::SharedPtr > const& readers
				) {
					this->process_one( changed, variants, readers, frequencyOutput, correlationOutput ) ;
				}
			) ;
			++count
		) {
			progress_context( count, visitor.count() ) ;
		} ;
	}
	
	void process_one(
		std::vector< int > const& changed,
		std::vector< genfile::VariantIdentifyingData > const& variants,
		std::vector< genfile::VariantDataReader::SharedPtr > const& readers,
		qcdb::Storage& frequencyOutput,
		qcdb::MultiVariantStorage& correlationOutput
	) {
		if( m_dosages.size() != changed.size() ) {
			m_dosages.resize( changed.size() ) ;
			m_ploidy.resize( changed.size() ) ;
			m_nonmissingness.resize( changed.size() ) ;
		}
		double const prior_weight = options().get< double >( "-prior-weight" ) ;
		for( std::size_t i = 0; i < changed.size(); ++i ) {
			if( changed[i] == 1 ) {
				DosageSetter setter( &(m_dosages[i]), &(m_ploidy[i]), &(m_nonmissingness[i]), options().check( "-assume-haploid" )) ;
				readers[i]->get( ":genotypes:", genfile::to_GP_unphased( setter )) ;

#if DEBUG
				std::cerr << globals::program_name + ":process_one(): variant " << i << ": " << variants[i] << ":\n" ;
				std::cerr << "Loaded data with nonmissingness: " << m_nonmissingness[i] << ".\n" ;
#endif

				FrequencyStore::const_iterator where = m_frequencies.find( variants[i] ) ;
				if( where == m_frequencies.end() ) {
					m_frequencies[ variants[i] ] = compute_regularised_frequency( m_dosages[i], m_ploidy[i], m_nonmissingness[i], prior_weight ) ;
					frequencyOutput.store_per_variant_data(
						variants[i],
						"frequency",
						m_frequencies[ variants[i] ]
					) ;
				}
			}
		}
		std::vector< metro::SampleRange > included_samples( 1, metro::SampleRange( 0, m_dosages[0].size() )) ;
		for( std::size_t i = 0; i < changed.size(); ++i ) {
			included_samples = metro::impl::intersect_ranges( included_samples, m_nonmissingness[i] ) ;
		}
		double covariance = 0.0, correlation = 0.0, N = 0.0 ;
		
		compute_regularised_correlation(
			m_dosages, m_ploidy, included_samples,
			&covariance, &correlation, &N,
			prior_weight
		) ;
		if( ( correlation * correlation ) >= options().get< double >( "-minimum-r2" ) ) {
			correlationOutput.store_data_for_key(
				variants,
				"N",
				int64_t(N)
			) ;
			// We store correlations as integers on the scale -16384 to 16384.
			// This is because sqlite has compression for integer storage, such that
			// reals take 8 bytes but integers take
			correlationOutput.store_data_for_key(
				variants,
				"encoded_r",
				// ints in sqlite use only 2 bytes if up to +ve integer 2287.
				// We map -1...1 to 0...2048, such that the transformation
				//
				// correlation = (stored_value / 1024) - 1.0
				// or
				// correlation = (stored_value - 1024.0) / 1024.0
				// maps back to correlation space.
				//
				int64_t( std::round((correlation + 1.0) * 1024 ))
			) ;
		}
#if DEBUG
		std::cerr << globals::program_name + ":process_one():\n" ;
		std::cerr << "dosage1: " << m_dosages[0].head( 15 ).transpose() << "...\n" ; 
		std::cerr << "dosage2: " << m_dosages[0].head( 15 ).transpose() << "...\n" ; 
		std::cerr << "N: " << N << ".\n" ;
		std::cerr << "covariance: " << covariance << ".\n" ;
		std::cerr << "correlation: " << correlation << ".\n" ;
		std::cerr << "included: " << included_samples << ".\n" ;
#endif
			
	}
	
private:
	std::vector< Eigen::VectorXd > m_dosages ;
	std::vector< Eigen::VectorXd > m_ploidy ;
	std::vector< std::vector< metro::SampleRange > > m_nonmissingness ;
	typedef std::map< genfile::VariantIdentifyingData, double > FrequencyStore ;
	FrequencyStore m_frequencies ;
	
private:
	
	double compute_regularised_frequency(
		Eigen::VectorXd const& dosages,
		Eigen::VectorXd const& ploidy,
		std::vector< metro::SampleRange > const& included_samples,
		double prior_weight = 1.0
	) {
		double result = 0.0 ;
		double N = 0.0 ;
		for( std::size_t i = 0; i < included_samples.size(); ++i ) {
			metro::SampleRange const& range = included_samples[i] ;
			result += dosages.segment( range.begin(), range.size() ).sum() ;
			N += ploidy.segment( range.begin(), range.size() ).sum() ;
		}
		// Assume we have additional data adding up to prior_weight haploid samples
		// evenly split between both alleles
		result += prior_weight / 2 ;
		N += prior_weight ;
		return result / N ;
	}

	void compute_regularised_correlation(
		std::vector< Eigen::VectorXd > const& dosages,
		std::vector< Eigen::VectorXd > const& ploidy,
		std::vector< metro::SampleRange > const& included_samples,
		double* covariance,
		double* correlation,
		double* number_of_samples,
		double const prior_weight = 1.0
	) {
		assert( dosages.size() == 2 ) ;
		assert( covariance ) ;
		assert( correlation ) ;
		assert( number_of_samples ) ;
		assert( ploidy[0] == ploidy[1] ) ;
		std::vector< double > frequencies( dosages.size(), 0.0 ) ;

		typedef Eigen::VectorBlock< Eigen::VectorXd const > ConstBlock ;
		typedef Eigen::VectorBlock< Eigen::VectorXi const > ConstIntegerBlock ;

		// We compute correlation augmented by additional data adding up to prior_weight
		// haploid samples evenly split between the four possible genotype combinations
		for( std::size_t v = 0; v < frequencies.size(); ++v ) {
			frequencies[v] = compute_regularised_frequency(
				dosages[v], ploidy[v], included_samples, prior_weight
			) ;
		}

		double result = 0.0 ;
		double N = 0.0 ;
		std::vector< double > variances( dosages.size(), 0.0 ) ;
		for( std::size_t i = 0; i < included_samples.size(); ++i ) {
			metro::SampleRange const& range = included_samples[i] ;
			ConstBlock d0 = dosages[0].segment( range.begin(), range.size() ) ;
			ConstBlock d1 = dosages[1].segment( range.begin(), range.size() ) ;
			ConstBlock p = ploidy[0].segment( range.begin(), range.size() ) ;
			result += (
				( d0 - p * frequencies[0] ).array()
				* ( d1 - p * frequencies[1] ).array()
			).sum() ;
			N += p.array().sum() ;
			for( std::size_t v = 0; v < frequencies.size(); ++v ) {
				ConstBlock d = dosages[v].segment( range.begin(), range.size() ) ;
				ConstBlock p = ploidy[v].segment( range.begin(), range.size() ) ;
				variances[v] += ( d - p * frequencies[v] ).array().square().sum() ;
			}
		}
		
		// Add regularising information consisting of a total of
		// prior_weight haploid samples, split evenly between
		// the four possible genotype combinations
		result += (0 - frequencies[0]) * (0 - frequencies[1] ) * (prior_weight / 4) ;
		result += (0 - frequencies[0]) * (1 - frequencies[1] ) * (prior_weight / 4) ;
		result += (1 - frequencies[0]) * (0 - frequencies[1] ) * (prior_weight / 4) ;
		result += (1 - frequencies[0]) * (1 - frequencies[1] ) * (prior_weight / 4) ;
		for( std::size_t v = 0; v < frequencies.size(); ++v ) {
			variances[v] += ( 0 - frequencies[v] ) * ( 0 - frequencies[v] ) * (prior_weight / 2.0) ;
			variances[v] += ( 1 - frequencies[v] ) * ( 1 - frequencies[v] ) * (prior_weight / 2.0) ;
		}
		N += prior_weight ;

		*covariance = result / N ;
		*correlation = result / (std::sqrt( variances[0] ) * std::sqrt( variances[1] )) ;
		*number_of_samples = N ;
	}
	
	struct DosageSetter: public genfile::VariantDataReader::PerSampleSetter {
		typedef std::vector< metro::SampleRange > SampleRanges ;

		DosageSetter(
			Eigen::VectorXd* dosages,
			Eigen::VectorXd* ploidy,
			SampleRanges* nonmissing_samples,
			bool treat_as_haploid = false
		):
			m_dosages( dosages ),
			m_ploidy( ploidy ),
			m_total_prob(0.0),
			m_treat_as_haploid( treat_as_haploid ),
			m_number_of_entries( 0 ),
			m_nonmissing_samples( nonmissing_samples ),
			m_sample_i(0)
		{
			assert( dosages != 0 && nonmissing_samples != 0 && ploidy != 0 ) ;
		}

		void initialise( std::size_t nSamples, std::size_t nAlleles ) {
			if( nAlleles != 2 ) {
				throw genfile::BadArgumentError(
					"ProbSetter::initialise()",
					"nAlleles=" + genfile::string_utils::to_string( nAlleles ),
					"I only support biallelic variants"
				) ;
			}
			m_dosages->resize( nSamples ) ;
			m_dosages->setZero() ;
			m_ploidy->resize( nSamples ) ;
			m_ploidy->setZero() ;
			
			m_nonmissing_samples->clear() ;
			m_last_nonmissing_sample_i = 0 ;
		}

		bool set_sample( std::size_t n ) {
			m_sample_i = n ;
			return true ;
		}

		void set_number_of_entries(
			uint32_t ploidy, std::size_t n,
			genfile::OrderType const order_type,
			genfile::ValueType const value_type
		) {
			if( value_type != genfile::eProbability ) {
				throw genfile::BadArgumentError(
					"Callsetter::set_number_of_entries()",
					"value_type=" + genfile::string_utils::to_string( value_type ),
					"Expected genotype call probabilities (e.g. GP field)."
				) ;
			}
			if( order_type != genfile::ePerUnorderedGenotype ) {
				throw genfile::BadArgumentError(
					"Callsetter::set_number_of_entries()",
					"order_type=" + genfile::string_utils::to_string( order_type ),
					"Expected genotype call probabilities (e.g. GP field)."
				) ;
			}
			if( m_treat_as_haploid ) {
				(*m_ploidy)(m_sample_i) = 1 ;
			} else {
				(*m_ploidy)(m_sample_i) = ploidy ;
			}
			m_number_of_entries = n ;
			m_total_prob = 0.0 ;
		}

		void set_value( std::size_t entry_i, genfile::MissingValue const value ) {
			// This sample has missing data, end our current sample range here.
			if( m_sample_i > m_last_nonmissing_sample_i ) {
				m_nonmissing_samples->push_back(
					metro::SampleRange( m_last_nonmissing_sample_i, m_sample_i )
				) ;
#if DEBUG > 2
				std::cerr << "DosageSetter::set_value(): added sample range: " << m_nonmissing_samples->back() << ".\n" ;
#endif
			}
			// Skip this sample for next range.
			m_last_nonmissing_sample_i = m_sample_i + 1 ;
		}

		void set_value( std::size_t entry_i, Integer const value ) {
			assert(0) ; // expected a probabilty
		}

		void set_value( std::size_t entry_i, double const value ) {
			if( m_sample_i >= m_last_nonmissing_sample_i ) {
				m_total_prob += value ;
				// handle missing data encoded as zeroes
				if( m_total_prob == 0.0 && (entry_i + 1) == m_number_of_entries ) {
					// missing
					m_nonmissing_samples->push_back(
						metro::SampleRange( m_last_nonmissing_sample_i, m_sample_i )
					) ;
					// Skip this sample for next range.
					m_last_nonmissing_sample_i = m_sample_i + 1 ;
				} else {
					// we assume two alleles, so entry_i is also the count of the B allele 
					// in the genotype.
					if( m_treat_as_haploid ) {
						// only homozygous genotypes are counted, others are assumed missing
						// entry_i == 0 contributes 0 dosage anyway
						if( entry_i == (m_number_of_entries - 1) ) {
							(*m_dosages)(m_sample_i) += value ;
						}
					} else {
						(*m_dosages)(m_sample_i) += value * entry_i ;
						m_total_prob += value ;
					}
				}
#if DEBUG > 2
				std::cerr << "DosageSetter::set_value(): added value " << (value * entry_i) << " to entry " << entry_i << " for sample " << m_sample_i << ".\n" ;
#endif
			}
		}

		void finalise() {
			// Add the final sample range if needed.
			if( (m_sample_i+1) > m_last_nonmissing_sample_i ) {
				m_nonmissing_samples->push_back(
					metro::SampleRange( m_last_nonmissing_sample_i, m_sample_i+1 )
				) ;
#if DEBUG > 2
				std::cerr << "DosageSetter::finalise(): added sample range: " << m_nonmissing_samples->back() << ".\n" ;
#endif
			}
#if DEBUG > 2
				std::cerr << "DosageSetter::finalise(): complete.\n" ;
#endif
		}
		
	private:
		Eigen::VectorXd* m_dosages ;
		Eigen::VectorXd* m_ploidy ;
		double m_total_prob ;
		bool m_treat_as_haploid ;

		std::size_t m_number_of_entries ;
		SampleRanges* m_nonmissing_samples ;
		std::size_t m_last_nonmissing_sample_i ;
		std::size_t m_sample_i ;
	} ;
} ;

int main( int argc, char** argv ) {
    try {
		CorrelatorApplication app( argc, argv ) ;
    }
	catch( appcontext::HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}

