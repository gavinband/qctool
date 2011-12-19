#include <string>
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread.hpp>
#include <Eigen/Core>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/RandomAccessSNPDataSource.hpp"
#include "genfile/SNPDataSourceChain.hpp"
#include "genfile/SampleFilteringSNPDataSource.hpp"
#include "genfile/wildcard.hpp"
#include "genfile/vcf/get_set.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "HaplotypeFrequencyComponent.hpp"
#include "integration/NewtonRaphson.hpp"
#include "integration/Derivative.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "db/Error.hpp"

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
		.set_default_value( "ld.sqlite3" ) ;
	options.option_implies_option( "-compute-ld-file", "-compute-ld-with" ) ;
	options.option_implies_option( "-compute-ld-with", "-cohort-name" ) ;
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
			std::string const& cohort,
			genfile::SNPIdentifyingData const& source_snp,
			genfile::SNPIdentifyingData const& target_snp,
			std::string const& variable,
			genfile::VariantEntry const& value
		) {
			if( !m_sink.get() ) {
				m_sink = statfile::BuiltInTypeStatSink::open( m_filename ) ;
				(*m_sink) | "cohort" | "SNPID1" | "rsid1" | "chromosome1" | "position1" | "SNPID2" | "variable" | "value" ;
			}
			
			(*m_sink)
				<< cohort
				<< source_snp.get_SNPID() << source_snp.get_rsid() << source_snp.get_position().chromosome() << source_snp.get_position().position()
				<< target_snp.get_SNPID() << target_snp.get_rsid() << target_snp.get_position().chromosome() << target_snp.get_position().position()
				<< variable << value.as< double >()
				<< statfile::end_row() ;
		}
		
	private:
		std::string const m_filename ;
		statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
	} ;
	
	struct HaplotypeFrequencyDBOutputter {
		typedef std::auto_ptr< HaplotypeFrequencyDBOutputter > UniquePtr ;
		typedef boost::shared_ptr< HaplotypeFrequencyDBOutputter > SharedPtr ;

		static UniquePtr create( std::string const& filename ) { return UniquePtr( new HaplotypeFrequencyDBOutputter( filename ) ) ; }
		static SharedPtr create_shared( std::string const& filename ) { return SharedPtr( new HaplotypeFrequencyDBOutputter( filename ) ) ; }

		HaplotypeFrequencyDBOutputter( std::string const& filename ):
			m_connection( db::Connection::create( filename )),
			m_max_transaction_count( 10000 )
		{
			db::Connection::ScopedTransactionPtr transaction = m_connection->open_transaction() ;
			m_connection->run_statement(
				"CREATE TABLE IF NOT EXISTS Variant ( id INTEGER PRIMARY KEY, snpid TEXT, rsid TEXT, chromosome TEXT, position INTEGER, alleleA TEXT, alleleB TEXT )"
			) ;
			m_connection->run_statement(
				"CREATE INDEX IF NOT EXISTS Variant_index ON Variant( chromosome, position, rsid )"
			) ;
			m_connection->run_statement(
				"CREATE TABLE IF NOT EXISTS Entity ( id INTEGER PRIMARY KEY, name TEXT, description TEXT )"
			) ;
			m_connection->run_statement(
				"CREATE TABLE IF NOT EXISTS PairwiseSummaryData ( "
				"cohort_id INT, variant1_id INT, variant2_id INT, variable_id INT, value FLOAT, "
				"FOREIGN KEY( cohort_id ) REFERENCES Entity( id ), "
				"FOREIGN KEY( variant1_id ) REFERENCES Variant( id ), "
				"FOREIGN KEY( variant2_id ) REFERENCES Variant( id ), "
				"FOREIGN KEY( variable_id ) REFERENCES Entity( id ), "
				"UNIQUE( cohort_id, variant1_id, variant2_id, variable_id ) "
				")"
			) ;
			m_connection->run_statement(
				"CREATE VIEW IF NOT EXISTS LDView AS "
				"SELECT C.name AS cohort, "
				"S1.chromosome AS chromosome1, S1.position AS position1, S1.rsid AS rsid1, "
				"S2.chromosome AS chromosome1, S2.position AS position1, S2.rsid AS rsid1, "
				"V.name AS variable, PSD.value AS value "
				"FROM PairwiseSummaryData PSD "
				"INNER JOIN Entity C ON C.id == PSD.cohort_id "
				"INNER JOIN Entity V ON V.id == PSD.variable_id "
				"INNER JOIN Variant S1 ON S1.id == PSD.variant1_id "
				"INNER JOIN Variant S2 ON S2.id == PSD.variant2_id "
			) ;
				
			m_connection->run_statement(
				"CREATE INDEX IF NOT EXISTS PairwiseSummaryDataIndex ON PairwiseSummaryData( cohort_id, variant1_id, variant2_id, variable_id )"
			) ;

			construct_statements() ;
		}

		~HaplotypeFrequencyDBOutputter() {
			write_data( m_data ) ;
		}

		void operator()(
			std::string const& cohort,
			genfile::SNPIdentifyingData const& source_snp,
			genfile::SNPIdentifyingData const& target_snp,
			std::string const& variable,
			genfile::VariantEntry const& value
		) {
			m_data.resize( m_data.size() + 1 ) ;
			m_data.back().get<0>() = cohort ;
			m_data.back().get<1>() = source_snp ;
			m_data.back().get<2>() = target_snp ;
			m_data.back().get<3>() = variable ;
			m_data.back().get<4>() = value ;

			if( m_data.size() == m_max_transaction_count ) {
				write_data( m_data ) ;
				m_data.clear() ;
			}
		}

	private:
		db::Connection::UniquePtr m_connection ;
		std::size_t const m_max_transaction_count ;

		db::Connection::StatementPtr m_find_variant_statement ;
		db::Connection::StatementPtr m_insert_variant_statement ;
		db::Connection::StatementPtr m_find_entity_statement ;
		db::Connection::StatementPtr m_insert_entity_statement ;
		db::Connection::StatementPtr m_insert_summarydata_statement ;

		typedef std::vector< boost::tuple< std::string, genfile::SNPIdentifyingData, genfile::SNPIdentifyingData, std::string, genfile::VariantEntry > > Data ;
		Data m_data ;

	private:
		void construct_statements() {
			m_find_variant_statement = m_connection->get_statement(
				"SELECT id FROM Variant WHERE rsid == ?1 AND chromosome == ?2 AND position == ?3"
			) ;
			m_insert_variant_statement = m_connection->get_statement(
				"INSERT INTO Variant ( snpid, rsid, chromosome, position, alleleA, alleleB) "
				"VALUES( ?1, ?2, ?3, ?4, ?5, ?6 )"
			) ;
			m_find_entity_statement = m_connection->get_statement( "SELECT * FROM Entity WHERE name == ?1" ) ;
			m_insert_entity_statement = m_connection->get_statement( "INSERT INTO Entity ( name ) VALUES ( ?1 )" ) ;
			m_insert_summarydata_statement = m_connection->get_statement(
				"INSERT OR REPLACE INTO PairwiseSummaryData ( cohort_id, variant1_id, variant2_id, variable_id, value ) "
				"VALUES( ?1, ?2, ?3, ?4, ?5 )"
			) ;
		}

		void write_data( Data const& data ) {
			db::Connection::ScopedTransactionPtr transaction ;

			for( std::size_t i = 0; i < 100; ++i ) {
				try {
					transaction = m_connection->open_transaction() ;
					break ;
				}
				catch( db::StatementStepError const& e ) {
					// wait a tenth of a second
					std::cerr << "HaplotypeFrequencyDBOutputter::write_data(): failed to open transaction, trying again in 0.1s...\n" ;
					boost::this_thread::sleep( boost::posix_time::milliseconds( 100 ) ) ;
				}
				catch( ... ) {
					std::cerr << "HaplotypeFrequencyDBOutputter::write_data(): OMG, a strange exception was caught.\n" ;
					boost::this_thread::sleep( boost::posix_time::milliseconds( 100 ) ) ;
				}
			}
			if( !transaction.get() ) {
				throw genfile::OperationFailedError( "HaplotypeFrequencyDBOutputter::write_data()", m_connection->get_spec(), "Opening transaction." ) ;
			}
			for( std::size_t i = 0; i < m_data.size(); ++i ) {
				store_comparison(
					data[i].get<0>(),
					data[i].get<1>(),
					data[i].get<2>(),
					data[i].get<3>(),
					data[i].get<4>()
				) ;
			}
		}

		db::Connection::RowId get_or_create_snp( genfile::SNPIdentifyingData const& snp ) const {
			m_find_variant_statement->reset()
				.bind( 1, snp.get_rsid() )
				.bind( 2, std::string( snp.get_position().chromosome() ))
				.bind( 3, snp.get_position().position() )
				.step()
			;
			if( m_find_variant_statement->empty() ) {
				m_insert_variant_statement
					->reset()
					.bind( 1, snp.get_SNPID() )
					.bind( 2, snp.get_rsid() )
					.bind( 3, std::string( snp.get_position().chromosome() ) )
					.bind( 4, snp.get_position().position() )
					.bind( 5, snp.get_first_allele())
					.bind( 6, snp.get_second_allele())
					.step()
				;

				return m_connection->get_last_insert_row_id() ;
			} else {
				return m_find_variant_statement->get< db::Connection::RowId >( 0 ) ;
			}
		}

		db::Connection::RowId get_or_create_entity( std::string const& name ) const {
			m_find_entity_statement
				->reset()
				.bind( 1, name ).step() ;

			if( m_find_entity_statement->empty() ) {
				m_insert_entity_statement
					->reset()
					.bind( 1, name )
					.step() ;
				return m_connection->get_last_insert_row_id() ;
			} else {
				return m_find_entity_statement->get< db::Connection::RowId >( 0 ) ;
			}
		}

		void store_comparison(
			std::string const& cohort,
			genfile::SNPIdentifyingData const& snp1,
			genfile::SNPIdentifyingData const& snp2,
			std::string const& variable,
			genfile::VariantEntry const& value
		) {

			db::Connection::RowId snp1_id = get_or_create_snp( snp1 ) ;
			db::Connection::RowId snp2_id = get_or_create_snp( snp2 ) ;
			db::Connection::RowId cohort_id = get_or_create_entity( cohort ) ;
			db::Connection::RowId variable_id = get_or_create_entity( variable );

			m_insert_summarydata_statement
				->reset()
				.bind( 1, cohort_id )
				.bind( 2, snp1_id )
				.bind( 3, snp2_id )
				.bind( 4, variable_id )
				.bind( 5, value.as< double >()  )
				.step()
			;
		}
	} ;
}

HaplotypeFrequencyComponent::UniquePtr HaplotypeFrequencyComponent::create(
	appcontext::OptionProcessor const& options,
	appcontext::UIContext& ui_context,
	std::vector< std::size_t > const& indices_of_filtered_out_samples
) {
	HaplotypeFrequencyComponent::UniquePtr result ;

	genfile::SNPDataSource::UniquePtr source = genfile::SNPDataSource::create_chain(
		genfile::wildcard::find_files_by_chromosome(
			options.get< std::string >( "-compute-ld-with" ),
			genfile::wildcard::eALL_CHROMOSOMES
		)
	) ;
	source = genfile::SampleFilteringSNPDataSource::create(
		source,
		std::set< std::size_t >( indices_of_filtered_out_samples.begin(), indices_of_filtered_out_samples.end() )
	) ;
	
	std::cerr << "LD SNP samples: " << source->number_of_samples() << " (filtered from " << source->get_parent_source().number_of_samples() << ").\n" ;
	
	result.reset(
		new HaplotypeFrequencyComponent(
			source,
			ui_context
		)
	) ;
	
	impl::HaplotypeFrequencyDBOutputter::SharedPtr outputter = impl::HaplotypeFrequencyDBOutputter::create_shared(
		options.get_value< std::string >( "-compute-ld-file" )
	) ;

	result->send_results_to(
		boost::bind(
			&impl::HaplotypeFrequencyDBOutputter::operator(),
			outputter,
			options.get< std::string >( "-cohort-name" ),
			_1,
			_2,
			_3,
			_4
		)
	) ;
	/*
	impl::HaplotypeFrequencyFileOutputter::SharedPtr outputter = impl::HaplotypeFrequencyFileOutputter::create_shared(
		options.get_value< std::string >( "-compute-ld-file" )
	) ;
	
	result->send_results_to(
		boost::bind(
			&impl::HaplotypeFrequencyFileOutputter::operator(),
			outputter,
			options.get< std::string >( "-cohort-name" ),
			_1,
			_2,
			_3,
			_4
		)
	) ;
	*/
	return result ;
}

HaplotypeFrequencyComponent::HaplotypeFrequencyComponent(
	genfile::SNPDataSource::UniquePtr source,
	appcontext::UIContext& ui_context
):
	m_source( source ),
	m_ui_context( ui_context ),
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
	try {
		compute_ld_measures(
			source_snp,
			target_snp,
			genotypes
		) ;
	} catch( genfile::OperationFailedError const& e ) {
		m_ui_context.logger() << "!! HaplotypeFrequencyComponent::compute_ld_measures(): could not compute LD measures between SNPs "
			<< source_snp << " and " << target_snp << ".\n"
			<< "!! reason: " << e.get_message() << "\n"
			<< "!! this comparison will not appear in output.\n" ;
	}
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
	HaplotypeFrequencyLogLikelihood::Vector pi( 4 ) ;
	pi.tail( 3 ) = ll.get_MLE_by_EM() ;
	pi(0) = 1.0 - pi.tail( 3 ).sum() ;
	double D = pi(0) * pi(3) - pi(1) * pi(2) ;
	double max_D ;
	if( D < 0 ) {
		max_D = std::min( (pi(0) + pi(1)) * (pi(0)+pi(2)), (pi(1)+pi(3))*(pi(2)+pi(3)) ) ;
	}
	else {
		max_D = std::min( (pi(0) + pi(1)) * (pi(1)+pi(3)), (pi(0)+pi(2))*(pi(2)+pi(3)) ) ;
	}
	double Dprime = D / max_D ;
	double r = D / std::sqrt( (pi(0)+pi(1)) * (pi(2)+pi(3)) * (pi(0)+pi(2)) * (pi(1)+pi(3))) ;

	m_result_signal( source_snp, target_snp, "pi00", pi(0) ) ;
	m_result_signal( source_snp, target_snp, "pi01", pi(1) ) ;
	m_result_signal( source_snp, target_snp, "pi10", pi(2) ) ;
	m_result_signal( source_snp, target_snp, "pi11", pi(3) ) ;
	m_result_signal( source_snp, target_snp, "D", D ) ;
	m_result_signal( source_snp, target_snp, "Dprime", Dprime ) ;
	m_result_signal( source_snp, target_snp, "r", r ) ;
	m_result_signal( source_snp, target_snp, "r_squared", r * r ) ;
}

void HaplotypeFrequencyComponent::end_processing_snps() {}

void HaplotypeFrequencyComponent::send_results_to( ResultCallback callback ) {
	m_result_signal.connect( callback ) ;
}
