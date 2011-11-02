#include <string>
#include <boost/signals2/signal.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "PairwiseCallComparer.hpp"
#include "PairwiseCallComparerManager.hpp"
#include "CallComparerComponent.hpp"

namespace impl {
	struct CallComparerFileOutputter {
		typedef std::auto_ptr< CallComparerFileOutputter > UniquePtr ;
		typedef boost::shared_ptr< CallComparerFileOutputter > SharedPtr ;
		
		static UniquePtr create( std::string const& filename ) { return UniquePtr( new CallComparerFileOutputter( filename ) ) ; }

		CallComparerFileOutputter( std::string const& filename ):
			m_filename( filename ),
			m_sink( statfile::BuiltInTypeStatSink::open( filename ))
		{
			(*m_sink) | "SNPID" | "rsid" | "chromosome" | "position" | "alleleA" | "alleleB" | "callset_1" | "callset_2" | "comparison_method" | "comparison_variable" | "value" ;
		}

		void write_comparison(
			genfile::SNPIdentifyingData const& snp,
			std::string const& callset1,
			std::string const& callset2,
			std::string const& comparison_method,
			std::string const& comparison_variable,
			genfile::VariantEntry const& value
		) {
			(*m_sink)
				<< snp.get_SNPID()
				<< snp.get_rsid()
				<< std::string( snp.get_position().chromosome() )
				<< snp.get_position().position()
				<< std::string( 1, snp.get_first_allele() )
				<< std::string( 1, snp.get_second_allele() )
				<< callset1
				<< callset2
				<< comparison_method
				<< comparison_variable ;
			(*m_sink)
				<< value.as< double >()
				<< statfile::end_row() ;
			;
		}

	private:
		std::string const m_filename ;
		statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
	} ;

	struct CallComparerDBOutputter {
		typedef std::auto_ptr< CallComparerDBOutputter > UniquePtr ;
		typedef boost::shared_ptr< CallComparerDBOutputter > SharedPtr ;
		
		static UniquePtr create( std::string const& filename ) { return UniquePtr( new CallComparerDBOutputter( filename ) ) ; }

		CallComparerDBOutputter( std::string const& filename ):
			m_connection( db::Connection::create( filename )),
			m_max_transaction_count( 10000 ),
			m_transaction_count( 0 )
		{
			m_connection->run_statement(
				"CREATE TABLE IF NOT EXISTS Variant ( id INTEGER PRIMARY KEY, snpid TEXT, rsid TEXT, chromosome TEXT, position INTEGER, alleleA TEXT, alleleB TEXT )"
			) ;
			m_connection->run_statement(
				"CREATE INDEX IF NOT EXISTS Variant_index ON Variant( rsid, chromosome, position )"
			) ;
			m_connection->run_statement(
				"CREATE TABLE IF NOT EXISTS Entity ( entity_id INTEGER PRIMARY KEY, name TEXT, description TEXT )"
			) ;
			m_connection->run_statement(
				"CREATE TABLE IF NOT EXISTS Comparison ( "
				"variant_id INT, callset1 TEXT, callset2 TEXT, method_id INT, variable_id INT, value FLOAT, "
				"FOREIGN KEY( variant_id ) REFERENCES Variant( id )), "
				"FOREIGN KEY( method_id ) REFERENCES Entity( id )), "
				"FOREIGN KEY( variable_id ) REFERENCES Entity( id ))"
			) ;
			m_connection->run_statement( "DELETE FROM Variant" ) ;
			m_connection->run_statement( "DELETE FROM Comparison" ) ;
			m_connection->run_statement( "DELETE FROM Entity" ) ;
			m_connection->run_statement( "BEGIN TRANSACTION" ) ;
		}
		
		~CallComparerDBOutputter() {
			m_connection->run_statement( "COMMIT TRANSACTION" ) ;
		}

		void write_comparison(
			genfile::SNPIdentifyingData const& snp,
			std::string const& callset1,
			std::string const& callset2,
			std::string const& comparison_method,
			std::string const& comparison_variable,
			genfile::VariantEntry const& value
		) {
			db::Connection::StatementPtr statement = m_connection->get_statement(
				"SELECT id FROM Variant WHERE rsid == ?1 AND chromosome == ?2 AND position == ?3"
			) ;
			statement
				->bind( 1, snp.get_rsid() )
				.bind( 2, std::string( snp.get_position().chromosome() ))
				.bind( 3, snp.get_position().position() )
				.step()
			;

			db::Connection::RowId snp_id ;

			if( statement->empty() ) {
				m_connection->get_statement(
					"INSERT INTO Variant ( snpid, rsid, chromosome, position, alleleA, alleleB) "
					"VALUES( ?1, ?2, ?3, ?4, ?5, ?6 )"
				)
					->bind( 1, snp.get_SNPID() )
					.bind( 2, snp.get_rsid() )
					.bind( 3, std::string( snp.get_position().chromosome() ) )
					.bind( 4, snp.get_position().position() )
					.bind( 5, std::string( 1, snp.get_first_allele() ))
					.bind( 6, std::string( 1, snp.get_second_allele() ))
					.step()
				;
				
				snp_id = m_connection->get_last_insert_row_id() ;
			} else {
				snp_id = statement->get< db::Connection::RowId >( 0 ) ;
			}

			db::Connection::RowId method_id ;
			db::Connection::RowId variable_id ;

			{
				statement = m_connection->get_statement( "SELECT * FROM Entity WHERE name == ?1" ) ;
				statement->bind( 1, comparison_method ).step() ;

				if( statement->empty() ) {
					m_connection->get_statement( "INSERT INTO Entity ( name ) VALUES ( ?1 )" )
						->bind( 1, comparison_method )
						.step() ;
					method_id = m_connection->get_last_insert_row_id() ;
				} else {
					method_id = statement->get< db::Connection::RowId >( 0 ) ;
				}

				statement = m_connection->get_statement( "SELECT * FROM Entity WHERE name == ?1" ) ;
				statement->bind( 1, comparison_variable ).step() ;

				if( statement->empty() ) {
					m_connection->get_statement( "INSERT INTO Entity ( name ) VALUES ( ?1 )" )
						->bind( 1, comparison_variable )
						.step() ;
					variable_id = m_connection->get_last_insert_row_id() ;
				} else {
					variable_id = statement->get< db::Connection::RowId >( 0 ) ;
				}
			}
			
			m_connection->get_statement(
				"INSERT INTO Comparison ( variant_id, callset1, callset2, method_id, variable_id, value ) "
				"VALUES( ?1, ?2, ?3, ?4, ?5, ?6 )"
			)
				->bind( 1, snp_id )
				.bind( 2, callset1 )
				.bind( 3, callset2 )
				.bind( 4, method_id )
				.bind( 5, variable_id )
				.bind( 6, value.as< double >()  )
				.step()
			;

			if( ++m_transaction_count == m_max_transaction_count ) {
				m_connection->run_statement( "COMMIT TRANSACTION" ) ;
				m_connection->run_statement( "BEGIN TRANSACTION" ) ;
				m_transaction_count = 0 ;
			}
		}

	private:
		db::Connection::UniquePtr m_connection ;
		std::size_t const m_max_transaction_count ;
		std::size_t m_transaction_count ;
	} ;
}

void CallComparerComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Call comparison options" ) ;
	options[ "-compare-calls" ]
		.set_description( "Compare genotype calls from the given fields of a VCF or other file. "
		 	"The value should be a comma-separated list of genotype call fields in the data.")
		.set_takes_single_value() ;
	options[ "-compare-calls-file" ]
		.set_description( "Name of output file to put call comparisons in." )
		.set_takes_single_value()
		.set_default_value( "call-comparisons.csv" ) ;
}

CallComparerComponent::UniquePtr CallComparerComponent::create( PairwiseCallComparerManager::UniquePtr comparer, std::vector< std::string > const& call_fields  ) {
	UniquePtr result(
		new CallComparerComponent( comparer, call_fields )
	) ;
	return result ;
}

CallComparerComponent::UniquePtr CallComparerComponent::create( appcontext::OptionProcessor const& options ) {
	PairwiseCallComparerManager::UniquePtr comparer( PairwiseCallComparerManager::create().release() ) ;

	std::string filename = options.get< std::string >( "-compare-calls-file" ) ;
	std::string const sqlite_suffix = ".sqlite3" ;
	if( filename.size() > sqlite_suffix.size() && ( filename.compare( filename.size() - sqlite_suffix.size(), sqlite_suffix.size(), sqlite_suffix) == 0 )) {
		comparer->send_results_to(
			boost::bind(
				&impl::CallComparerDBOutputter::write_comparison,
				impl::CallComparerDBOutputter::SharedPtr( impl::CallComparerDBOutputter::create( filename )),
				_1, _2, _3, _4, _5, _6
			)
		) ;
	} else {
		comparer->send_results_to(
			boost::bind(
				&impl::CallComparerFileOutputter::write_comparison,
				impl::CallComparerFileOutputter::SharedPtr( impl::CallComparerFileOutputter::create( filename )),
				_1, _2, _3, _4, _5, _6
			)
		) ;
	}

	comparer->add_comparer(
		"AlleleFrequencyTestCallComparer",
		PairwiseCallComparer::create( "AlleleFrequencyTestCallComparer" )
	) ;

	return CallComparerComponent::create(
		comparer,
		genfile::string_utils::split_and_strip_discarding_empty_entries(
			options.get_value< std::string >( "-compare-calls" ),
			",",
			" \t"
		)
	) ;
}

CallComparerComponent::CallComparerComponent( PairwiseCallComparerManager::UniquePtr call_comparer, std::vector< std::string > const& call_fields  ):
	m_call_comparer( call_comparer ),
	m_call_fields( call_fields )
{}

void CallComparerComponent::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
}

void CallComparerComponent::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	// Add all the calls to the call comparer.
	m_call_comparer->set_SNP( snp ) ;
	genfile::SingleSNPGenotypeProbabilities calls ;
	for( std::size_t i = 0; i < m_call_fields.size(); ++i ) {
		data_reader.get( m_call_fields[i], calls ) ;
		m_call_comparer->add_calls( m_call_fields[i], calls ) ;
	}
}

void CallComparerComponent::end_processing_snps() {
}
