#include <iostream>
#include <iomanip>
#include "genfile/string_utils.hpp" 
#include "genfile/Error.hpp" 
#include "genfile/zlib.hpp" 
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "VCDBDataStore.hpp"

VCDBDataStore::VCDBDataStore( db::Connection::UniquePtr connection ):
	m_connection( connection ),
	m_zlib_compression_id( 0 ),
	m_no_compression_id( 0 )
{
	setup_db( *m_connection ) ;
	load_variants( *m_connection ) ;
}

void VCDBDataStore::setup_db( db::Connection& connection ) {
	db::Connection::StatementPtr statement ;
	std::string version ;
	try {
		// Let's see if the database is empty (a new database).
		statement = connection.get_statement(
			"SELECT * FROM sqlite_master"
		) ;
		statement->step() ;
		if( statement->empty() ) {
			m_db_version = "VCDBv0.3" ;
			prepare_new_db( connection, m_db_version ) ;
			std::cerr << "VCDBDataStore: Prepared new DB with format \"" << m_db_version << "\".\n" ;
		}
		else {
			m_db_version = get_db_version( connection ) ;
			prepare_existing_db( connection, m_db_version ) ;
			std::cerr << "VCDBDataStore: Opened DB with format " << m_db_version << ".\n" ;
		}
	}
	catch( db::Error const& e ) {
		std::cerr << "Database " << m_connection->get_spec() << " does not appear to be in vcdb format (no VCDBFileVersion table).\n" ;
		throw genfile::MalformedInputError( m_connection->get_spec(), 0 ) ;
	}
}

void VCDBDataStore::load_variants( db::Connection& connection ) {
	db::Connection::StatementPtr statement = connection.get_statement(
		"SELECT id, rsid, chromosome, position, alleleA, alleleB FROM SNP"
	) ;
	for( statement->step() ; !statement->empty(); statement->step() ) {
		m_variants.resize( m_variants.size() + 1 ) ;
		m_variants.back().first.rsid() = statement->get< std::string >(1) ;
		m_variants.back().first.SNPID() = statement->get< std::string >(1) ;
		m_variants.back().first.position().chromosome() = genfile::Chromosome( statement->get< std::string >(2) ) ;
		m_variants.back().first.position().position() = statement->get< int64_t >( 3 ) ;
		m_variants.back().first.first_allele() = statement->get< std::string >( 4 )[0] ;
		m_variants.back().first.second_allele() = statement->get< std::string >( 5 )[0] ;
	}
	std::sort( m_variants.begin(), m_variants.end() ) ;
	std::cerr << "Loaded " << m_variants.size() << " variants.\n" ;
}

std::string VCDBDataStore::get_db_version( db::Connection& connection ) const {
	db::Connection::StatementPtr statement ;
	statement = connection.get_statement(
		"SELECT key, value FROM FileInfo WHERE key == ?1"
	) ;
	statement->bind( 1, "format" ).step() ;
	if( statement->empty() ) {
		throw genfile::MalformedInputError( m_connection->get_spec() + ":FileInfo", 0 ) ;
	}
	std::string version = statement->get< std::string >( 1 ) ;
	if( version != "VCDBv0.2" && version != "VCDBv0.3" ) {
		std::cerr << "Database " << m_connection->get_spec() << ": version \"" + version + "\" is not supported.\n" ;
		throw genfile::MalformedInputError( m_connection->get_spec() + ":FileInfo", 0 ) ;
	}
	return version ;
}

void VCDBDataStore::prepare_existing_db( db::Connection& connection, std::string const& version ) {
	prepare_db( connection, version, false ) ;
}

void VCDBDataStore::prepare_new_db( db::Connection& connection, std::string const& version ) {
	prepare_db( connection, version, true ) ;
}

void VCDBDataStore::prepare_db( db::Connection& connection, std::string const& version, bool new_db ) {
	db::Connection::StatementPtr statement ;

	std::string stub = new_db ? "" : "IF NOT EXISTS " ;

	if( new_db ) {
		connection.run_statement(
			"PRAGMA page_size = 8192"
		) ;
		connection.run_statement(
			"CREATE TABLE FileInfo ( key TEXT NOT NULL UNIQUE, value TEXT );"
		) ;
		statement = connection.get_statement(
			"INSERT INTO FileInfo VALUES( ?1, ?2 )"
		) ;
		statement->bind( 1, "format" ).bind( 2, "VCDBv0.3" ).step() ;
	}

	connection.run_statement(
		"CREATE TABLE " + stub + "Entity ( id INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE, description TEXT )"
	) ;

	connection.run_statement(
		"CREATE TABLE " + stub + "EntityRelationship ( "
			"entity1_id INTEGER NOT NULL, "
			"relationship_id INTEGER NOT NULL, "
			"entity2_id INTEGER NOT NULL, "
			"UNIQUE ( entity1_id, entity2_id, relationship_id ), "
			"FOREIGN KEY( entity1_id ) REFERENCES Entity( id ), "
			"FOREIGN KEY( entity2_id ) REFERENCES Entity( id ), "
			"FOREIGN KEY( relationship_id ) REFERENCES Entity( id ) "
		");"
	) ;

	connection.run_statement(
		"CREATE VIEW " + stub + "EntityRelationshipNameView AS "
		"SELECT ER.entity1_id, E1.name AS entity1_name, ER.relationship_id, ERel.name AS relationship_name, ER.entity2_id, E2.name AS entity2_name "
		"FROM EntityRelationship ER "
		"INNER JOIN Entity E1 ON ER.entity1_id = E1.id "
		"INNER JOIN Entity ERel ON ER.relationship_id = ERel.id "
		"INNER JOIN Entity E2 ON ER.entity2_id = E2.id ;"
	) ;

	connection.run_statement(
		"CREATE TABLE " + stub + "SNP ( "
		"id INTEGER PRIMARY KEY, "
		"rsid TEXT NOT NULL, "
		"chromosome TEXT NOT NULL, "
		"position INTEGER NOT NULL, "
		"alleleA TEXT NOT NULL, "
		"alleleB TEXT NOT NULL, "
		"UNIQUE( rsid, chromosome, position ) "
		")"
	) ;
	
	connection.get_statement(
		"CREATE TABLE " + stub + "VariantLocation ("
		"snp_id INTEGER NOT NULL, "
		"build_id INTEGER NOT NULL, "
		"chromosome TEXT NOT NULL, "
		"position INTEGER NOT NULL, "
		"alleleA TEXT NOT NULL, "
		"alleleB TEXT NOT NULL, "
		"FOREIGN KEY (snp_id) REFERENCES SNP( id ), "
		"FOREIGN KEY (snp_id) REFERENCES Entity( id )"
		")"
	) ;
	
	if( m_db_version >= "VCDBv0.3" ) {
		connection.run_statement(
			"CREATE TABLE " + stub + "Data ( "
				"snp_id INTEGER NOT NULL, "
				"field_id INTEGER NOT NULL, "
				"cohort_id INTEGER NOT NULL, "
				"storage_id INTEGER NOT NULL, "
				"compression_id INTEGER NOT NULL, "
				"uncompressed_size INTEGER NOT NULL, "
				"data BLOB, "
				"UNIQUE( snp_id, field_id, cohort_id, storage_id ), "
				"FOREIGN KEY (snp_id) REFERENCES SNP( id ), "
				"FOREIGN KEY (field_id) REFERENCES Entity( id ), "
				"FOREIGN KEY (cohort_id) REFERENCES Entity( id ), "
				"FOREIGN KEY (storage_id) REFERENCES Entity( id ), "
				"FOREIGN KEY (compression_id) REFERENCES Entity( id ) "
			") ;"
		) ;
	} else {
		connection.run_statement(
			"CREATE TABLE " + stub + "Data ( "
				"snp_id INTEGER NOT NULL, "
				"field_id INTEGER NOT NULL, "
				"cohort_id INTEGER NOT NULL, "
				"storage_id INTEGER NOT NULL, "
				"uncompressed_size INTEGER NOT NULL, "
				"data BLOB, "
				"UNIQUE( snp_id, field_id, cohort_id, storage_id ), "
				"FOREIGN KEY (snp_id) REFERENCES SNP( id ), "
				"FOREIGN KEY (field_id) REFERENCES Entity( id ), "
				"FOREIGN KEY (cohort_id) REFERENCES Entity( id ), "
				"FOREIGN KEY (storage_id) REFERENCES Entity( id )"
			") ;"
		) ;
	}

	connection.run_statement(
		"CREATE VIEW " + stub + "DataNameView AS "
		"SELECT S.*, D.field_id, F.id AS field_Id, F.name AS field_name, C.id AS cohort_id, C.name AS cohort_name, ST.id AS storage_id, ST.name AS storage_name, uncompressed_size FROM Data D "
		"INNER JOIN SNP S ON S.id == D.snp_id "
		"INNER JOIN Entity F ON F.id == D.field_id "
		"INNER JOIN Entity C ON C.id == D.cohort_id "
		"INNER JOIN Entity ST ON ST.id == D.storage_id "
	) ;
	
	get_or_create_entity( "is_a", "Indicates that a is an object of class b." ) ;
	get_or_create_entity( "stored_as", "Indicates that a is stored in the format specified by b." ) ;
	get_or_create_entity( "has_type", "Indicates that values of field a have type b." ) ;
	get_or_create_entity( "storage_type", "Class of entities representing storage types." ) ;
	if( m_db_version <= "VCDBv0.2" ) {
		get_or_create_entity( "zlib_compressed_serialized", "zlib compressed data, serialized using qctool." ) ;
		set_relationship( "zlib_compressed_serialized", "is_a", "storage_type" ) ;
	}
	else {
		get_or_create_entity( "compression_type", "Class of entities representing compression algorithms." ) ;
		m_zlib_compression_id = get_or_create_entity( "zlib", "zlib compressed data" ) ;
		m_no_compression_id = get_or_create_entity( "uncompressed", "uncompressed data" ) ;
		set_relationship( "zlib", "is_a", "compression_type" ) ;
		set_relationship( "uncompressed", "is_a", "compression_type" ) ;
	}
	get_or_create_entity( "cohort", "A cohort" ) ;
	get_or_create_entity( "per_variant_per_sample_data", "Indicates data pertains to each sample at a variant" ) ;
	get_or_create_entity( "per_variant_data", "Indicates data pertains to a variant as a whole" ) ;
	get_or_create_entity( "unnamed cohort", "An unnamed cohort" ) ;
	
	create_indices( connection ) ;
}

void VCDBDataStore::create_indices( db::Connection& connection ) {
	std::string const stub = "IF NOT EXISTS " ;
	connection.run_statement(
		"CREATE INDEX " + stub + "Entity_name ON Entity( name )"
	) ;

	connection.run_statement(
		"CREATE INDEX " + stub + "SNP_rsid ON SNP( rsid )"
	) ;
	connection.run_statement(
		"CREATE INDEX " + stub + "SNP_position ON SNP( chromosome, position )"
	) ;

	connection.run_statement(
		"CREATE INDEX " + stub + "Data_snp ON Data( snp_id )"
	) ;
	connection.run_statement(
		"CREATE INDEX " + stub + "Data_field ON Data( field_id )"
	) ;
	connection.run_statement(
		"CREATE INDEX " + stub + "Data_storage ON Data( storage_id )"
	) ;
	connection.run_statement(
		"CREATE INDEX " + stub + "Data_cohort ON Data( cohort_id )"
	) ;
}

std::string VCDBDataStore::get_spec() const {
	return "VCDBDataStore:" + m_connection->get_spec() ;
}

int64_t VCDBDataStore::get_or_create_SNP( genfile::SNPIdentifyingData const& snp ) {
	using genfile::string_utils::to_string ;
	
	db::Connection::StatementPtr statement  = m_connection->get_statement(
		"SELECT id, alleleA, alleleB FROM SNP WHERE rsid = ?1 AND chromosome = ?2 AND position = ?3"
	) ;
	statement
		->bind( 1, snp.get_rsid() )
		.bind( 2, to_string( snp.get_position().chromosome() ) )
		.bind( 3, snp.get_position().position() ) ;
	statement->step() ;

	db::Connection::RowId result ;

	if( statement->empty() ) {
		statement  = m_connection->get_statement(
			"INSERT INTO SNP( rsid, chromosome, position, alleleA, alleleB ) VALUES( ?1, ?2, ?3, ?4, ?5 )"
			) ;
		statement
			->bind( 1, snp.get_rsid() )
			.bind( 2, to_string( snp.get_position().chromosome() ))
			.bind( 3, snp.get_position().position() )
			.bind( 4, std::string( 1, snp.get_first_allele() ) )
			.bind( 5, std::string( 1, snp.get_second_allele() ) ) ;
		statement->step() ;
		result = m_connection->get_last_insert_row_id() ;
	}
	else {
		result = statement->get_column< db::Connection::RowId >( 0 ) ;
		std::string alleleA = statement->get_column< std::string >( 1 ) ;
		std::string alleleB = statement->get_column< std::string >( 2 ) ;
		if( alleleA != std::string( 1, snp.get_first_allele() ) || alleleB != std::string( 1, snp.get_second_allele() )) {
			throw genfile::MismatchError(
				"VCDBDataStore::get_or_create_SNP()", m_connection->get_spec() + ":SNP",
				to_string( snp ),
				"entry with id " + to_string( result )
			) ;
		}
	}
	return result ;
}

int64_t VCDBDataStore::get_or_create_entity( std::string name, std::string description ) {
	db::Connection::RowId result ;
	db::Connection::StatementPtr statement = m_connection->get_statement(
		"SELECT id, description FROM Entity WHERE name == ?1"
	) ;
	statement->bind( 1, name ) ;
	statement->step() ;
	if( statement->empty() ) {
		statement = m_connection->get_statement(
			"INSERT INTO Entity ( name, description ) VALUES ( ?1, ?2 )"
		) ;
		statement->bind( 1, name ) ;
		statement->bind( 2, description ) ;
		statement->step() ;
		result = m_connection->get_last_insert_row_id() ;
	} else if( !statement->is_null( 1 )) {
		std::string const stored_description = statement->get_column< std::string >( 1 ) ;
		if( stored_description != description ) {
			throw genfile::BadArgumentError( "IntensityWriter::get_or_create_entity()", "description = \"" + description + "\"" ) ;
		}
		result = statement->get_column< db::Connection::RowId >( 0 ) ;
	}
	return result ;
}

int64_t VCDBDataStore::get_entity( std::string const& name ) const {
	db::Connection::RowId result ;
	
	db::Connection::StatementPtr statement = m_connection->get_statement(
		"SELECT id, description FROM Entity WHERE name == ?1"
	) ;
	statement->bind( 1, name ) ;
	statement->step() ;
	if( statement->empty() ) {
		throw genfile::BadArgumentError( "VCDBDataStore::get_entity()", "VCDBDataStore::get_entity(): no entity with name \"" + name + "\" exists.\n" ) ;
	}
	else {
		result = statement->get< db::Connection::RowId >( 0 ) ;
	}
	return result ;
}

void VCDBDataStore::set_relationship( db::Connection::RowId left_id, db::Connection::RowId relation_id, db::Connection::RowId right_id ) const {
	db::Connection::StatementPtr statement = m_connection->get_statement(
		"INSERT OR REPLACE INTO EntityRelationship VALUES( ?1, ?2, ?3 )"
	) ;
	statement
		->bind( 1, left_id )
		.bind( 2, relation_id )
		.bind( 3, right_id )
	;
	statement->step() ;
}

void VCDBDataStore::set_relationship( std::string const& left, std::string const& relation, std::string const& right ) const {
	db::Connection::RowId
		left_id = get_entity( left ),
		relation_id = get_entity( relation ),
		right_id = get_entity( right ) ;
		
	set_relationship( left_id, relation_id, right_id ) ;
}

void VCDBDataStore::store_per_variant_data( int64_t snp_id, std::string const& field, std::string const& cohort, std::string const& storage, char const* buffer, char const* const end ) {
	db::Connection::RowId
		field_id = get_entity( field ),
		cohort_id = get_entity( cohort ),
		storage_id = get_entity( storage ) ;

	store_per_variant_data( snp_id, field_id, cohort_id, storage_id, buffer, end ) ;
}

db::Connection::StatementPtr VCDBDataStore::get_store_variant_data_statement( std::string const& db_version ) {
	std::string SQL = "INSERT OR REPLACE INTO Data (snp_id, field_id, cohort_id, storage_id, " ;
	if( db_version >= "VCDBv0.3" ) {
		SQL += "uncompressed_size, data, compression_id ) VALUES( ?1, ?2, ?3, ?4, ?5, ?6, ?7 )" ;
	} else {
		SQL += "uncompressed_size, data ) VALUES( ?1, ?2, ?3, ?4, ?5, ?6 )" ;
	}
	return m_connection->get_statement( SQL ) ;
}

void VCDBDataStore::store_per_variant_data(
	int64_t snp_id,
	int64_t field_id,
	int64_t cohort_id,
	int64_t storage_id,
	char const* buffer,
	char const* const end
) {
	assert( end >= buffer ) ;

	if( !m_store_variant_data_statement.get() ) {
		m_store_variant_data_statement = get_store_variant_data_statement( m_db_version ) ;
	}
	else {
		m_store_variant_data_statement->reset() ;
	}
	db::Connection::StatementPtr& statement = m_store_variant_data_statement ;
	statement->bind( 1, snp_id ) ;
	statement->bind( 2, field_id ) ;
	statement->bind( 3, cohort_id ) ;
	statement->bind( 4, storage_id ) ;
	statement->bind( 5, int64_t( end - buffer ) ) ;

	bool compress_data = (m_db_version <= "VCDBv0.2") || (( end - buffer ) > 1000) ;
	if( compress_data ) {
		// Compress the data before storing.
		genfile::zlib_compress( buffer, end, &m_compression_buffer ) ;
		statement->bind( 6, &m_compression_buffer[0], &m_compression_buffer[0] + m_compression_buffer.size() 	) ;
		/*
		std::cerr << "First few bytes of compressed data are:\n" ;
		for( std::size_t i = 0; i < 10; ++i ) {
			std::cerr << " " << std::hex << std::setw(2) << std::setfill( '0' ) << unsigned( static_cast< unsigned char >( m_compression_buffer[i] ) ) ;
		}
		std::cerr << std::dec << "...\n" ;
		*/
	}
	else {
		statement->bind( 6, buffer, end ) ;
	}
	if( m_db_version >= "VCDBv0.3" ) {
		statement->bind( 7, ( compress_data ? m_zlib_compression_id : m_no_compression_id ) ) ;
	}
	statement->step() ;
}

void VCDBDataStore::get_entities_by_relation(
	std::string const& relation,
	std::string const& related_entity,
	boost::function< void ( db::Connection::RowId, std::string const& ) > callback
) {
	get_entities_by_relation( get_entity( relation ), get_entity( related_entity ), callback ) ;
}

void VCDBDataStore::get_entities_by_relation(
	int64_t relation_id,
	int64_t related_entity_id,
	boost::function< void ( db::Connection::RowId, std::string const& ) > callback
) {
	db::Connection::StatementPtr statement = m_connection->get_statement(
		"SELECT E.id, E.name FROM Entity E INNER JOIN EntityRelationship ER ON ER.entity1_id == E.id AND ER.relation_id = ?1 AND ER.entity2_id = ?2"
	) ;
	statement
		->bind( 1, relation_id )
		.bind( 2, related_entity_id ) ;
	for( statement->step(); !statement->empty(); statement->step() ) {
		callback( statement->get< db::Connection::RowId >( 0 ), statement->get< std::string >( 1 )) ;
	}
}

DataStore::TransactionPtr VCDBDataStore::open_transaction() {
	return DataStore::TransactionPtr( new VCDBTransaction( *m_connection ) ) ;
}

VCDBDataStore::VCDBTransaction::VCDBTransaction( db::Connection& connection ):
	m_connection( &connection ),
	m_rollback( false )
{
	m_connection->run_statement( "BEGIN" ) ;
}

VCDBDataStore::VCDBTransaction::~VCDBTransaction() {
	if( m_rollback ) {
		m_connection->run_statement( "ROLLBACK" ) ;
	}
	else {
		m_connection->run_statement( "COMMIT" ) ;
	}
}

