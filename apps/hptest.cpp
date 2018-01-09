
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
#include "db/Error.hpp"
#include "metro/RegressionDesign.hpp"

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
		
		options.declare_group( "Model options" ) ;
		options[ "-model" ]
			.set_description( "Model to fit" )
			.set_takes_values(1)
			.set_default_value( "gp ~ gh" )
		;
	}
} ;


namespace {
	struct CallSetter: public genfile::VariantDataReader::PerSampleSetter {
		CallSetter(
			Eigen::VectorXd* genotypes,
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
			m_genotypes->resize( nSamples ) ;
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
				// accumulate count of the 2nd allele
				(*m_genotypes)(m_sample_i) += value ;
			}
		}

		void set_value( std::size_t entry_i, double const value ) {
			assert(0) ; // expecting GT field
		}

		void finalise() {
			// nothing to do
		}
		
	private:
		Eigen::VectorXd* m_genotypes ;
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
		
		test( *host, *para, *samples ) ;
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
		genfile::CohortIndividualSource const& samples
	) {
		// We model the predictor as a fixed covariate
		std::vector< std::string > covariate_names ;
		std::vector< std::string > predictor_names ;
		assert( options().get< std::string >( "-model" ) == "gp ~ gh" ) ;
		covariate_names.push_back( "gh" ) ;

		metro::RegressionDesign::UniquePtr design = metro::RegressionDesign::create(
			Vector::Zero( samples.size() ), Vector::Constant( samples.size(), 1.0 ), "gp",
			Matrix::Zero( samples.size(), 1 ), Matrix::Zero( samples.size(), 1 ), covariate_names,
			predictor_names
		) ;
		
		genfile::SNPDataSource* const predictor_source = &host ;
		genfile::SNPDataSource* const outcome_source = &para ;
		
		// Make the predictor the inner loop
		genfile::VariantIdentifyingData ov, pv ;
		Vector outcome = Vector::Zero( samples.size() ) ;
		Vector outcome_nonmissingness = Vector::Constant( samples.size(), 1.0 ) ;
		Vector outcome_ploidy = Vector::Zero( samples.size() ) ;
		Vector predictor = Vector::Zero( samples.size() ) ;
		Vector predictor_nonmissingness = Vector::Constant( samples.size(), 1.0 ) ;
		Vector predictor_ploidy = Vector::Zero( samples.size() ) ;
		while( predictor_source->get_snp_identifying_data( &pv )) {
			predictor_source->read_variant_data()->get( "GT", CallSetter( &predictor, &predictor_nonmissingness, &predictor_ploidy )) ;
			design->set_covariates( predictor, predictor_nonmissingness ) ;
			outcome_source->reset_to_start() ;
			while( outcome_source->get_snp_identifying_data( &ov )) {
				outcome_source->read_variant_data()->get( "GT", CallSetter( &outcome, &outcome_nonmissingness, &outcome_ploidy ) ) ;
				design->set_outcome( outcome, outcome_nonmissingness, "gp" ) ;
				
				get_ui_context().logger()
					<< "test(): read data for outcome variant: ["
					<< ov
					<< "] and predictor variant ["
					<< pv
					<< "]:\n"
					<< design->get_summary()
					<< "\n" ;
				;
			}
		}
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

