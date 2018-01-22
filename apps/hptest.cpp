
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <iostream>
#include "appcontext/CmdLineOptionProcessor.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "config/qctool_version_autogenerated.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPDataSourceChain.hpp"
#include "genfile/ToGP.hpp"
#include "db/Error.hpp"
#include "metro/regression/Design.hpp"
#include "metro/SampleRange.hpp"
#include "metro/regression/BinomialLogistic.hpp"
#include "metro/ValueStabilisesStoppingCondition.hpp"
#include "metro/Snptest25StoppingCondition.hpp"
#include "metro/maximisation.hpp"
#include "qcdb/MultiVariantStorage.hpp"
#include "qcdb/FlatTableDBOutputter.hpp"
#include <Eigen/Core>
#include <Eigen/QR>

namespace globals {
	std::string const program_name = "hptest" ;
	std::string const program_version = qctool_version ;
	std::string const program_revision =  std::string( qctool_revision ).substr( 0, 7 ) ;
}

struct HPTestOptions: public appcontext::CmdLineOptionProcessor
{
public:
	std::string get_program_name() const { return globals::program_name ; }

	void declare_options( appcontext::OptionProcessor& options ) {
		// Meta-options
		options.set_help_option( "-help" ) ;
		options.set_spec_option( "-spec" ) ;

		// File options
		options.declare_group( "Input file options" ) ;
	    options[ "-gh" ]
	        .set_description( 	"Path to host genotype files."
								"The given filename may contain the wildcard character '#', which expands to match a"
								"one- or two-character chromosome identifier." )
			.set_takes_values( 1 )
			.set_minimum_multiplicity( 1 )
			.set_maximum_multiplicity( 1 ) ;

		options[ "-gp" ]
			.set_description(
				"Path to parasite genotype files."
				"The given filename may contain the wildcard character '#', which expands to match a"
				"one- or two-character chromosome identifier." )
			.set_takes_values( 1 )
			.set_minimum_multiplicity( 1 )
			.set_maximum_multiplicity( 1 ) ;

	    options[ "-ghfiletype" ]
			.set_description(
				"Specify the filetype of the genotype files specified by -gh. "
				"By default, qctool will guess the file type.  Use this option to override that guess. "
				"Possible types are: \"" + genfile::string_utils::join( genfile::SNPDataSource::get_file_types(), "\",\"" ) + "\"." )
			.set_takes_single_value()
			.set_default_value( "guess" ) ;

	    options[ "-gpfiletype" ]
			.set_description(
				"Specify the filetype of the genotype files specified by -gp. "
				"By default, qctool will guess the file type.  Use this option to override that guess. "
				"Possible types are: \"" + genfile::string_utils::join( genfile::SNPDataSource::get_file_types(), "\",\"" ) + "\"." )
			.set_takes_single_value()
			.set_default_value( "guess" ) ;
		
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
		options[ "-model" ]
			.set_description( "Model to fit" )
			.set_takes_values(1)
			.set_default_value( "gp ~ gh" )
		;

		options[ "-tolerance" ]
			.set_description( "Tolerance" )
			.set_takes_values(1)
			.set_default_value( 0.001 )
		;

		options[ "-max-iterations" ]
			.set_description( "Maximum fitting iterations" )
			.set_takes_values(1)
			.set_default_value( 100 )
		;
		
		options.declare_group( "Misc options" ) ;
		
		options[ "-compare-variants-by" ]
			.set_description( "By default, matching SNPs between cohorts uses all the available fields"
				" (position, rsid, snpid, and alleles.)"
				" Use this option to specify a comma-separated subset of those fields to use instead."
				" The first entry must be \"position\"."
				" This option can be used, for example, when cohorts are typed on different platforms so have different SNPID fields." )
			.set_takes_single_value()
			.set_default_value( "position,alleles,ids" ) ;
		options[ "-analysis-name" ]
			.set_description( "Specify a name to label results from this analysis with." )
			.set_takes_single_value()
			.set_default_value( "qctool analysis" ) ;
		options[ "-analysis-chunk" ]
			.set_description( "Specify a name denoting the current genomic region or chunk on which this is run.  This is intended for use in parallel environments." )
			.set_takes_single_value()
			.set_default_value( genfile::MissingValue() ) ;
	}
} ;


namespace {
	struct CallSetter: public genfile::VariantDataReader::PerSampleSetter {
		CallSetter(
			Eigen::MatrixXd* genotypes,
			Eigen::VectorXd* nonmissingness,
			Eigen::VectorXd* ploidy
		):
			m_genotypes( genotypes ),
			m_nonmissingness( nonmissingness ),
			m_ploidy( ploidy ),
			m_sample_i(0)
		{
			assert( genotypes != 0 && nonmissingness != 0 && ploidy != 0 ) ;
		}

		void initialise( std::size_t nSamples, std::size_t nAlleles ) {
			if( nAlleles != 2 ) {
				throw genfile::BadArgumentError(
					"CallSetter::initialise()",
					"nAlleles=" + genfile::string_utils::to_string( nAlleles ),
					"I only support biallelic variants"
				) ;
			}
			m_genotypes->resize( nSamples, nAlleles ) ;
			m_nonmissingness->resize( nSamples ) ;
			m_ploidy->resize( nSamples ) ;
			
			m_genotypes->setZero() ;
			m_nonmissingness->setZero() ;
			m_ploidy->setZero() ;
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
			if( value_type != genfile::eAlleleIndex ) {
				throw genfile::BadArgumentError(
					"Callsetter::set_number_of_entries()",
					"value_type=" + genfile::string_utils::to_string( value_type ),
					"Expected a hard-called genotype, consider using threshholded calls."
				) ;
			}
			if( order_type != genfile::ePerOrderedHaplotype && order_type != genfile::ePerUnorderedHaplotype ) {
				throw genfile::BadArgumentError(
					"Callsetter::set_number_of_entries()",
					"order_type=" + genfile::string_utils::to_string( order_type ),
					"Expected a hard-called genotype, consider using threshholded calls."
				) ;
			}
			(*m_ploidy)(m_sample_i) = ploidy ;
			(*m_nonmissingness)(m_sample_i) = 1 ;
		}

		void set_value( std::size_t entry_i, genfile::MissingValue const value ) {
			// Nothing to do, value is already marked missing
			(*m_nonmissingness)(m_sample_i) = 0 ;
		}

		void set_value( std::size_t entry_i, Integer const value ) {
			// Only accumulate if not missing
			double const& nonmissing = (*m_nonmissingness)(m_sample_i) ;
			if( nonmissing == 1 ) {
				// accumulate count of either allele
				(*m_genotypes)(m_sample_i,value) += 1 ;
			}
		}

		void set_value( std::size_t entry_i, double const value ) {
			assert(0) ; // expecting GT field
		}

		void finalise() {
			// nothing to do
		}
		
	private:
		Eigen::MatrixXd* m_genotypes ;
		Eigen::VectorXd* m_nonmissingness ;
		Eigen::VectorXd* m_ploidy ;
		std::size_t m_sample_i ;
	} ;
	
	struct ProbSetter: public genfile::VariantDataReader::PerSampleSetter {
		ProbSetter(
			Eigen::MatrixXd* genotypes,
			Eigen::VectorXd* nonmissingness,
			Eigen::VectorXd* ploidy
		):
			m_genotypes( genotypes ),
			m_nonmissingness( nonmissingness ),
			m_ploidy( ploidy ),
			m_sample_i(0)
		{
			assert( genotypes != 0 && nonmissingness != 0 && ploidy != 0 ) ;
		}

		void initialise( std::size_t nSamples, std::size_t nAlleles ) {
			if( nAlleles != 2 ) {
				throw genfile::BadArgumentError(
					"CallSetter::initialise()",
					"nAlleles=" + genfile::string_utils::to_string( nAlleles ),
					"I only support biallelic variants"
				) ;
			}
			m_genotypes->resize( nSamples, 3 ) ;
			m_nonmissingness->resize( nSamples ) ;
			m_ploidy->resize( nSamples ) ;
			
			m_genotypes->setZero() ;
			m_nonmissingness->setZero() ;
			m_ploidy->setZero() ;
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
			(*m_nonmissingness)(m_sample_i) = 1 ;
		}

		void set_value( std::size_t entry_i, genfile::MissingValue const value ) {
			// Nothing to do, value is already marked missing
			(*m_nonmissingness)(m_sample_i) = 0 ;
			m_genotypes->row(m_sample_i).setZero() ;
		}

		void set_value( std::size_t entry_i, Integer const value ) {
			assert(0) ; // expected a probabilty
		}

		void set_value( std::size_t entry_i, double const value ) {
			double const& nonmissing = (*m_nonmissingness)(m_sample_i) ;
			if( nonmissing == 1 ) {
				(*m_genotypes)(m_sample_i, entry_i) = value ;
			}
		}

		void finalise() {
			// nothing to do
		}
		
	private:
		Eigen::MatrixXd* m_genotypes ;
		Eigen::VectorXd* m_nonmissingness ;
		Eigen::VectorXd* m_ploidy ;
		std::size_t m_sample_i ;
	} ;
}

struct HPTestApplication: public appcontext::ApplicationContext
{
public:
	
	typedef Eigen::MatrixXd Matrix ;
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::RowVectorXd RowVector ;
	
public:
	HPTestApplication( int argc, char** argv ):
		appcontext::ApplicationContext(
			globals::program_name,
			globals::program_version + ", revision " + globals::program_revision,
			std::auto_ptr< appcontext::OptionProcessor >( new HPTestOptions ),
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
		catch( db::Error const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << ") with the following statement: \""
				<< e.sql()
				<< "\".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
	}
	
	void unsafe_process() {
		using genfile::string_utils::to_string ;
		
		genfile::SNPDataSource::UniquePtr
			host = open_snp_data_sources( options().get< std::string >( "-gh" )) ;
		genfile::SNPDataSource::UniquePtr
			para = open_snp_data_sources( options().get< std::string >( "-gp" )) ;
		genfile::CohortIndividualSource::UniquePtr
			samples = genfile::CohortIndividualSource::create( options().get< std::string >( "-s" ) ) ;

		std::size_t const N = samples->size() ;
		if( host->number_of_samples() != N ) {
			throw genfile::BadArgumentError(
				"HPTestApplication::unsafe_process()",
				"-gh \"" + options().get< std::string >( "-gh" ) + "\"",
				"Wrong number of samples (" + to_string(host->number_of_samples()) + ", expected " + to_string(N) + ")"
			) ;
		}
		if( para->number_of_samples() != N ) {
			throw genfile::BadArgumentError(
				"HPTestApplication::unsafe_process()",
				"-gp \"" + options().get< std::string >( "-gp" ) + "\"",
				"Wrong number of samples (" + to_string(para->number_of_samples()) + ", expected " + to_string(N) + ")"
			) ;
		}
		
		write_preamble( *host, *para, *samples ) ;
		
		{
			qcdb::MultiVariantStorage::UniquePtr storage = qcdb::MultiVariantStorage::create(
				options().get< std::string > ( "-o" ),
				2,
				options().get< std::string > ( "-analysis-name" ),
				options().get< std::string > ( "-analysis-chunk" ),
				options().get_values_as_map()
			) ;
			storage->set_variant_names( std::vector< std::string >({ "predictor", "outcome" })) ;
			test( *host, *para, *samples, *storage ) ;
			storage->finalise() ;
		}
	}
	
	genfile::SNPDataSource::UniquePtr open_snp_data_sources( std::string const& filename ) {
		std::vector< genfile::wildcard::FilenameMatch > filenames
			= genfile::wildcard::find_files_by_chromosome(
				filename,
				genfile::wildcard::eALL_CHROMOSOMES
			) ;

		genfile::SNPDataSourceChain::UniquePtr source = genfile::SNPDataSourceChain::create(
			filenames
		) ;
		
		return genfile::SNPDataSource::UniquePtr( source.release() ) ;
	}

	void write_preamble(
		genfile::SNPDataSource const& host,
		genfile::SNPDataSource const& para,
		genfile::CohortIndividualSource const& samples
	) {
		get_ui_context().logger() << "Loaded data for " << samples.size() << " samples.\n" ;
		get_ui_context().logger() << "    Host data:\n" << host.get_summary() << "\n" ;
		get_ui_context().logger() << "Parasite data:\n" << para.get_summary() << "\n" ;
		get_ui_context().logger() << "\n" << "Fitting model: " << options().get< std::string >( "-model" ) ;
	}
	
	void test(
		genfile::SNPDataSource& host,
		genfile::SNPDataSource& para,
		genfile::CohortIndividualSource const& samples,
		qcdb::MultiVariantStorage& output
	) {
		typedef metro::regression::Design::Matrix Matrix ;
		using namespace metro ;

		// Put together data needed for setting up regression design
		assert( options().get< std::string >( "-model" ) == "gp ~ gh" ) ;

		std::vector< std::string > const outcome_names({ "gp=1", "gp=0" }) ;
		Matrix outcome = Matrix::Zero( samples.size(), 2 ) ;
		Vector outcome_nonmissingness = Vector::Constant( samples.size(), 1.0 ) ;
		Vector outcome_ploidy = Vector::Zero( samples.size() ) ;

		std::vector< std::string > covariate_names({}) ;
		Matrix covariates = Matrix::Zero( samples.size(), 0 ) ;
		Vector covariate_nonmissingness = Vector::Constant( samples.size(), 1.0 ) ;

		std::vector< std::string > predictor_names({ "gh" }) ;
		Matrix predictor_probabilities = Matrix::Zero( samples.size(), 2 ) ;
		Vector predictor_nonmissingness = Vector::Constant( samples.size(), 1.0 ) ;
		Vector predictor_ploidy = Vector::Zero( samples.size() ) ;
		Matrix predictor_levels = Matrix::Zero( 3, 1 ) ;
		predictor_levels.col(0) = Eigen::VectorXd::LinSpaced( 3, 0, 2 ) ;

		regression::Design::UniquePtr design = regression::Design::create(
			outcome, outcome_nonmissingness, outcome_names,
			covariates, covariate_nonmissingness, covariate_names,
			predictor_names
		) ;
		
		genfile::SNPDataSource* const predictor_source = &host ;
		genfile::SNPDataSource* const outcome_source = &para ;

		// Reserve some space for our variants
		std::vector< genfile::VariantIdentifyingData > variants( 2 ) ;
		genfile::VariantIdentifyingData& pv = variants[0] ;
		genfile::VariantIdentifyingData& ov = variants[1] ;

		while( predictor_source->get_snp_identifying_data( &pv )) {
			ProbSetter hostSetter( &predictor_probabilities, &predictor_nonmissingness, &predictor_ploidy ) ;
			predictor_source->read_variant_data()->get( "GP", genfile::to_GP_unphased( hostSetter )) ;
			design->set_predictors( predictor_levels, predictor_probabilities, design->nonmissing_samples() ) ;

			outcome_source->reset_to_start() ;
			get_ui_context().logger()
				<< "Predictors:\n"
				<< predictor_probabilities.block(0,0, std::min( 10, int( predictor_probabilities.rows() )), predictor_probabilities.cols() )
				<< "\n" ;

			outcome_source->reset_to_start() ;
			while( outcome_source->get_snp_identifying_data( &ov )) {
				outcome_source->read_variant_data()->get( "GT", CallSetter( &outcome, &outcome_nonmissingness, &outcome_ploidy ) ) ;

				get_ui_context().logger()
					<< "Outcome:\n"
					<< outcome.block(0,0, std::min( 10, int( outcome.rows() )), outcome.cols() )
					<< "\n" ;

				design->set_outcome( outcome, outcome_nonmissingness, std::vector< std::string >({"gp=0","gp=1"})) ;
				
				get_ui_context().logger()
					<< "test(): read data for outcome variant: ["
					<< ov
					<< "] and predictor variant ["
					<< pv
					<< "]:\n"
					<< design->get_summary()
					<< "\n" ;
				
				regression::BinomialLogistic::UniquePtr ll = regression::BinomialLogistic::create( *design ) ;
				get_ui_context().logger()
					<< "test(): loglikelihood is:\n"
						<< ll->get_summary() << "\n" ;
				;

				std::cerr << "Included samples = " << design->nonmissing_samples() << ".\n" ;
				
				Eigen::VectorXd const parameters = fitit( *ll ) ;
				
				{
					double const fpv
						= ((2 * predictor_probabilities.col(2).sum() ) + predictor_probabilities.col(1).sum())
						/ (2 * predictor_probabilities.sum() ) ;
					double const fov = outcome.col(1).sum() / outcome.sum() ;
			
					get_ui_context().logger()
						<< "###\n"
						<< "++ predictor variant: " << pv << ", frequency " << fpv << "\n"
						<< "++   outcome variant: " << ov << ", frequency " << fov << "\n"
						<< "++        parameters: " << parameters.transpose() << ".\n" ;
					
					output.create_new_key( variants ) ;
					output.store_data_for_key(
						variants, "predictor_frequency", fpv
					) ;
					output.store_data_for_key(
						variants, "outcome_frequency", fov
					) ;
				}
				
			}
		}
	}
	
	Eigen::VectorXd fitit( metro::regression::LogLikelihood& ll ) {
		// Let's optimise...
		using namespace metro ;
		Snptest25StoppingCondition< regression::LogLikelihood > stopping_condition(
			ll,
			options().get< double >( "-tolerance" ),
			options().get< std::size_t >( "-max-iterations" ),
			&get_ui_context().logger()
		) ;
		typedef metro::regression::Design::Matrix Matrix ;
		Eigen::ColPivHouseholderQR< Matrix > solver ;
		//Eigen::VectorXd parameters = Eigen::VectorXd::Zero(2) ;
		Eigen::VectorXd parameters = Eigen::VectorXd::Zero(2) ;
		return maximise_by_newton_raphson( ll, parameters, stopping_condition, solver ) ;
	}
} ;

int main( int argc, char** argv ) {
    try {
		HPTestApplication app( argc, argv ) ;
    }
	catch( appcontext::HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}

