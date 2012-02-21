#include "genfile/Error.hpp"
#include "components/CallComparerComponent/CallComparerComponent.hpp"
#include "components/CallComparerComponent/CallComparerDBOutputter.hpp"

CallComparerDBOutputter::UniquePtr CallComparerDBOutputter::create( std::string const& filename ) {
	return UniquePtr( new CallComparerDBOutputter( filename ) ) ;
}

CallComparerDBOutputter::SharedPtr CallComparerDBOutputter::create_shared( std::string const& filename ) {
	return SharedPtr( new CallComparerDBOutputter( filename ) ) ;
}

CallComparerDBOutputter::CallComparerDBOutputter( std::string const& filename ):
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
		"CREATE TABLE IF NOT EXISTS Comparison ( "
		"variant_id INT, callset1 INT, callset2 INT, method_id INT, variable_id INT, value FLOAT, "
		"FOREIGN KEY( variant_id ) REFERENCES Variant( id ), "
		"FOREIGN KEY( method_id ) REFERENCES Entity( id ), "
		"FOREIGN KEY( variable_id ) REFERENCES Entity( id ))"
	) ;
	m_connection->run_statement(
		"CREATE INDEX IF NOT EXISTS ComparisonIndex ON Comparison( variant_id, method_id, variable_id )"
	) ;

	m_connection->run_statement(
		"CREATE TABLE IF NOT EXISTS SummaryData ( "
		"variant_id INT, analysis_id INT, variable_id INT, value FLOAT, "
		"FOREIGN KEY( variant_id ) REFERENCES Variant( id ), "
		"FOREIGN KEY( analysis_id ) REFERENCES Entity( id ), "
		"FOREIGN KEY( variable_id ) REFERENCES Entity( id ), "
		"UNIQUE( analysis_id, variant_id, variable_id ) "
		")"
	) ;

	m_connection->run_statement(
		"CREATE INDEX IF NOT EXISTS SummaryDataIndex ON SummaryData( analysis_id, variant_id, variable_id )"
	) ;

	m_connection->run_statement(
		"CREATE VIEW IF NOT EXISTS SummaryDataView AS "
		"SELECT          V.id AS variant_id, V.chromosome, V.position, V.rsid, "
		"SD.analysis_id, Analysis.name, Variable.id AS variable_id, Variable.name AS variable, "
		"SD.value AS value "
		"FROM SummaryData SD "
		"INNER JOIN Variant V ON V.id == SD.variant_id "
		"INNER JOIN Entity Analysis ON Analysis.id = SD.analysis_id "
		"INNER JOIN Entity Variable ON Variable.id = SD.variable_id"
	) ;

	construct_statements() ;
}

CallComparerDBOutputter::~CallComparerDBOutputter() {
	write_data( m_data ) ;
}

void CallComparerDBOutputter::begin_comparisons( genfile::SNPIdentifyingData const& snp ) {
	m_snp = snp ;
}

void CallComparerDBOutputter::end_comparisons() {}

void CallComparerDBOutputter::set_result(
	std::string const& callset1,
	std::string const& callset2,
	std::string const& comparison_method,
	std::string const& comparison_variable,
	genfile::VariantEntry const& value
) {
	m_data.resize( m_data.size() + 1 ) ;
	m_data.back().get<0>() = m_snp ;
	m_data.back().get<1>() = callset1 ;
	m_data.back().get<2>() = callset2 ;
	m_data.back().get<3>() = comparison_method ;
	m_data.back().get<4>() = comparison_variable ;
	m_data.back().get<5>() = value ;

	if( m_data.size() == m_max_transaction_count ) {
		write_data( m_data ) ;
		m_data.clear() ;
	}
}

void CallComparerDBOutputter::set_result(
	std::string const& comparison_method,
	std::string const& comparison_variable,
	genfile::VariantEntry const& value
) {
	m_data.resize( m_data.size() + 1 ) ;
	m_data.back().get<0>() = m_snp ;
	m_data.back().get<1>() = "" ;
	m_data.back().get<2>() = "" ;
	m_data.back().get<3>() = comparison_method ;
	m_data.back().get<4>() = comparison_variable ;
	m_data.back().get<5>() = value ;

	if( m_data.size() == m_max_transaction_count ) {
		write_data( m_data ) ;
		m_data.clear() ;
	}
}

void CallComparerDBOutputter::construct_statements() {
	m_find_variant_statement = m_connection->get_statement(
		"SELECT id FROM Variant WHERE rsid == ?1 AND chromosome == ?2 AND position == ?3"
	) ;
	m_insert_variant_statement = m_connection->get_statement(
		"INSERT INTO Variant ( snpid, rsid, chromosome, position, alleleA, alleleB) "
		"VALUES( ?1, ?2, ?3, ?4, ?5, ?6 )"
	) ;
	m_find_entity_statement = m_connection->get_statement( "SELECT * FROM Entity WHERE name == ?1" ) ;
	m_insert_entity_statement = m_connection->get_statement( "INSERT INTO Entity ( name ) VALUES ( ?1 )" ) ;
	m_insert_comparison_statement = m_connection->get_statement(
		"INSERT INTO Comparison ( variant_id, callset1, callset2, method_id, variable_id, value ) "
		"VALUES( ?1, ?2, ?3, ?4, ?5, ?6 )"
	) ;
	m_insert_summarydata_statement = m_connection->get_statement(
		"INSERT OR REPLACE INTO SummaryData ( variant_id, analysis_id, variable_id, value ) "
		"VALUES( ?1, ?2, ?3, ?4 )"
	) ;
}

void CallComparerDBOutputter::write_data( Data const& data ) {
	db::Connection::ScopedTransactionPtr transaction ;

	for( std::size_t i = 0; i < 100; ++i ) {
		try {
			transaction = m_connection->open_transaction() ;
			break ;
		}
		catch( db::StatementStepError const& e ) {
			// wait a tenth of a second
			std::cerr << "CallComparerComponent::write_data(): failed to open transaction, trying again in 0.1s...\n" ;
			boost::this_thread::sleep( boost::posix_time::milliseconds( 100 ) ) ;
		}
		catch( ... ) {
			std::cerr << "CallComparerComponent::write_data(): OMG, a strange exception was caught.\n" ;
			boost::this_thread::sleep( boost::posix_time::milliseconds( 100 ) ) ;
		}
	}
	if( !transaction.get() ) {
		throw genfile::OperationFailedError( "CallComparerComponent::write_data()", m_connection->get_spec(), "Opening transaction." ) ;
	}
	for( std::size_t i = 0; i < m_data.size(); ++i ) {
		store_comparison(
			data[i].get<0>(),
			data[i].get<1>(),
			data[i].get<2>(),
			data[i].get<3>(),
			data[i].get<4>(),
			data[i].get<5>()
		) ;
	}
}

void CallComparerDBOutputter::store_comparison(
	genfile::SNPIdentifyingData const& snp,
	std::string const& callset1,
	std::string const& callset2,
	std::string const& comparison_method,
	std::string const& comparison_variable,
	genfile::VariantEntry const& value
) {
	m_find_variant_statement
		->bind( 1, snp.get_rsid() )
		.bind( 2, std::string( snp.get_position().chromosome() ))
		.bind( 3, snp.get_position().position() )
		.step()
	;

	db::Connection::RowId snp_id ;

	if( m_find_variant_statement->empty() ) {
		m_insert_variant_statement
			->bind( 1, snp.get_SNPID() )
			.bind( 2, snp.get_rsid() )
			.bind( 3, std::string( snp.get_position().chromosome() ) )
			.bind( 4, snp.get_position().position() )
			.bind( 5, snp.get_first_allele())
			.bind( 6, snp.get_second_allele())
			.step()
		;
		
		snp_id = m_connection->get_last_insert_row_id() ;
	} else {
		snp_id = m_find_variant_statement->get< db::Connection::RowId >( 0 ) ;
	}
	
	m_find_variant_statement->reset() ;
	m_insert_variant_statement->reset() ;

	db::Connection::RowId method_id ;
	db::Connection::RowId callset1_id ;
	db::Connection::RowId callset2_id ;
	db::Connection::RowId variable_id ;

	if( callset1 == "" ) {
		assert( callset2 == "" ) ;
		m_find_entity_statement
			->bind( 1, callset1 ).step() ;

		if( m_find_entity_statement->empty() ) {
			m_insert_entity_statement
				->bind( 1, callset1 )
				.step() ;
			callset1_id = m_connection->get_last_insert_row_id() ;
		} else {
			callset1_id = m_find_entity_statement->get< db::Connection::RowId >( 0 ) ;
		}

		m_find_entity_statement->reset()
			.bind( 1, callset2 ).step() ;

		if( m_find_entity_statement->empty() ) {
			m_insert_entity_statement->reset()
				.bind( 1, callset2 )
				.step() ;
			callset2_id = m_connection->get_last_insert_row_id() ;
		} else {
			callset2_id = m_find_entity_statement->get< db::Connection::RowId >( 0 ) ;
		}
	}

	m_find_entity_statement->reset()
		.bind( 1, comparison_method ).step() ;

	if( m_find_entity_statement->empty() ) {
		m_insert_entity_statement->reset()
			.bind( 1, comparison_method )
			.step() ;
		method_id = m_connection->get_last_insert_row_id() ;
	} else {
		method_id = m_find_entity_statement->get< db::Connection::RowId >( 0 ) ;
	}

	m_find_entity_statement->reset()
		.bind( 1, comparison_variable ).step() ;

	if( m_find_entity_statement->empty() ) {
		m_insert_entity_statement->reset()
			.bind( 1, comparison_variable )
			.step() ;
		variable_id = m_connection->get_last_insert_row_id() ;
	} else {
		variable_id = m_find_entity_statement->get< db::Connection::RowId >( 0 ) ;
	}
	
	m_find_entity_statement->reset() ;
	m_insert_entity_statement->reset() ;
	
	if( callset1 != "" ) {
		m_insert_comparison_statement
			->bind( 1, snp_id )
			.bind( 2, callset1_id )
			.bind( 3, callset2_id )
			.bind( 4, method_id )
			.bind( 5, variable_id )
			.bind( 6, value  )
			.step()
		;
		m_insert_comparison_statement->reset() ;
		
	} else {
		m_insert_summarydata_statement
			->bind( 1, snp_id )
			.bind( 2, method_id )
			.bind( 3, variable_id )
			.bind( 4, value  )
			.step() ;
		m_insert_summarydata_statement->reset() ;
	}
	
}
