#include <iostream>
#include <boost/bind.hpp>
#include "genfile/FileUtils.hpp"
#include "db/SQLite3Connection.hpp"
#include "db/SQLStatement.hpp"
#include "genfile/Error.hpp"
#include "genfile/zlib.hpp"
#include "genfile/endianness_utils.hpp"
#include "IntensityWriter.hpp"

IntensityWriter::IntensityWriter( std::string const& filename ):
	m_filename( filename ),
	m_connection( new db::SQLite3Connection( filename ))
{
	setup() ;
}

void IntensityWriter::setup() {
	db::Connection& connection = *m_connection ;
	db::Connection::StatementPtr statement ;
	connection.run_statement(
		"CREATE TABLE IF NOT EXISTS FileInfo ( key TEXT NOT NULL UNIQUE, value TEXT );"
	) ;
	statement = connection.get_statement(
		"INSERT OR REPLACE INTO FileInfo VALUES( ?1, ?2 )"
	) ;
	statement->bind( 1, "format" ) ;
	statement->bind( 2, "VCDBv0.2" ) ;
	statement->step() ;

	connection.run_statement(
		"CREATE TABLE IF NOT EXISTS Entity ( id INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE, description TEXT ) ; "
		"CREATE INDEX IF NOT EXISTS Entity_name ON Entity( name )"
	) ;

	get_or_create_entity( "is_a", "Indicates that a is an object of class b." ) ;
	get_or_create_entity( "stored_as", "Indicates that a is stored in the format specified by b." ) ;
	get_or_create_entity( "has_type", "Indicates that values of field a have type b." ) ;
	get_or_create_entity( "storage_type", "Class of entities representing storage types." ) ;
	get_or_create_entity( "zlib_compressed_serialized", "zlib compressed data, serialized using qctool." ) ;
	get_or_create_entity( "double", "A single floating-point literal in native format" ) ;
	get_or_create_entity( "cohort", "A cohort" ) ;
	get_or_create_entity( "per_variant_per_sample_data", "Indicates data pertains to each sample at a variant" ) ;
	get_or_create_entity( "per_variant_data", "Indicates data pertains to a variant as a whole" ) ;
	get_or_create_entity( "unnamed cohort", "An unnamed cohort" ) ;

	connection.run_statement(
		"CREATE TABLE IF NOT EXISTS EntityRelationship ( "
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
		"CREATE VIEW IF NOT EXISTS EntityRelationshipNameView AS "
		"SELECT ER.entity1_id, E1.name AS entity1_name, ER.relationship_id, ERel.name AS relationship_name, ER.entity2_id, E2.name AS entity2_name "
		"FROM EntityRelationship ER "
		"INNER JOIN Entity E1 ON ER.entity1_id = E1.id "
		"INNER JOIN Entity ERel ON ER.relationship_id = ERel.id "
		"INNER JOIN Entity E2 ON ER.entity2_id = E2.id ;"
	) ;

	set_relationship( "zlib_compressed_serialized", "is_a", "storage_type" ) ;
	set_relationship( "double", "is_a", "storage_type" ) ;

	connection.run_statement(
		"CREATE TABLE IF NOT EXISTS SNP ( "
		"id INTEGER PRIMARY KEY, "
		"rsid TEXT NOT NULL, "
		"chromosome TEXT NOT NULL, "
		"position INTEGER NOT NULL, "
		"alleleA TEXT NOT NULL, "
		"alleleB TEXT NOT NULL, "
		"UNIQUE( rsid, chromosome, position ) "
		")"
	) ;

	connection.run_statement(
		"CREATE INDEX IF NOT EXISTS SNP_rsid ON SNP( rsid )"
	) ;
	statement = connection.get_statement(
		"CREATE INDEX IF NOT EXISTS SNP_position ON SNP( chromosome, position )"
	) ;
	statement->step() ;
	
	statement = connection.get_statement(
		"CREATE TABLE IF NOT EXISTS Data ( "
			"snp_id INTEGER NOT NULL, "
			"field_id INTEGER NOT NULL, "
			"cohort_id INTEGER NOT NULL, "
			"storage_id INTEGER NOT NULL, "
			"uncompressed_size INTEGER NOT NULL, "
			"data BLOB, "
			"UNIQUE( snp_id, field_id, storage_id ), "
			"FOREIGN KEY (snp_id) REFERENCES SNP( id ), "
			"FOREIGN KEY (field_id) REFERENCES Entity( id ), "
			"FOREIGN KEY (cohort_id) REFERENCES Entity( id ), "
			"FOREIGN KEY (storage_id) REFERENCES Entity( id ) "
		") ;"
	) ;
	statement->step() ;
	connection.run_statement(
		"CREATE INDEX IF NOT EXISTS Data_snp ON Data( snp_id )"
	) ;
	connection.run_statement(
		"CREATE INDEX IF NOT EXISTS Data_field ON Data( field_id )"
	) ;
	connection.run_statement(
		"CREATE INDEX IF NOT EXISTS Data_storage ON Data( storage_id )"
	) ;
	connection.run_statement(
		"CREATE INDEX IF NOT EXISTS Data_cohort ON Data( cohort_id )"
	) ;
}

void IntensityWriter::get_or_create_entity( std::string const& name, std::string const& description ) {
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
	m_entities[ name ] = result ;
}

void IntensityWriter::set_relationship( std::string const& left, std::string const& relation, std::string const& right ) const {
	db::Connection::StatementPtr statement = m_connection->get_statement(
		"INSERT OR REPLACE INTO EntityRelationship VALUES( ?1, ?2, ?3 )"
	) ;
	std::map< std::string, db::Connection::RowId >::const_iterator where = m_entities.find( left ) ;
	assert( where != m_entities.end() ) ;
	statement->bind( 1, where->second ) ;
	where = m_entities.find( relation ) ;
	assert( where != m_entities.end() ) ;
	statement->bind( 2, where->second ) ;
	where = m_entities.find( right ) ;
	assert( where != m_entities.end() ) ;
	statement->bind( 3, where->second ) ;
	statement->step() ;
}


void IntensityWriter::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	m_number_of_samples = number_of_samples ;
	m_number_of_snps = number_of_snps ;
	m_number_of_snps_written = 0 ;
}

namespace impl {
	void add_spec_to_map( std::map< std::string, std::string >* map, std::string const& name, std::string const& type ) {
		(*map)[ name ] = type ;
	}
}

void IntensityWriter::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	std::map< std::string, std::string > fields ;
	data_reader.get_supported_specs( boost::bind( &impl::add_spec_to_map, &fields, _1, _2 )) ;
	
	try {
		db::Connection::StatementPtr statement ;
		db::Connection::StatementPtr transaction = m_connection->get_statement(
			"BEGIN"
		) ;
		transaction->step() ;

		db::Connection::RowId snp_row_id ;
		{
			statement  = m_connection->get_statement(
				"SELECT id, rsid, chromosome, position, alleleA, alleleB FROM SNP WHERE rsid = ?1 AND chromosome = ?2 AND position = ?3"
			) ;
			statement->bind( 1, snp.get_rsid() ) ;
			statement->bind( 2, genfile::string_utils::to_string( snp.get_position().chromosome() )) ;
			statement->bind( 3, snp.get_position().position() ) ;
			statement->step() ;
			if( statement->empty() ) {
				statement  = m_connection->get_statement(
					"INSERT INTO SNP( rsid, chromosome, position, alleleA, alleleB ) VALUES( ?1, ?2, ?3, ?4, ?5 )"
				) ;
				statement->bind( 1, snp.get_rsid() ) ;
				statement->bind( 2, genfile::string_utils::to_string( snp.get_position().chromosome() )) ;
				statement->bind( 3, snp.get_position().position() ) ;
				statement->bind( 4, std::string( 1, snp.get_first_allele() ) ) ;
				statement->bind( 5, std::string( 1, snp.get_second_allele() ) ) ;
				statement->step() ;
				snp_row_id = m_connection->get_last_insert_row_id() ;
			}
			else {
				snp_row_id = statement->get_column< db::Connection::RowId >( 0 ) ;
				std::string alleleA = statement->get_column< std::string >( 4 ) ;
				std::string alleleB = statement->get_column< std::string >( 5 ) ;
				if( alleleA != std::string( 1, snp.get_first_allele() ) || alleleB != std::string( 1, snp.get_second_allele() )) {
					throw genfile::MismatchError(
						"IntensityWriter::processed_snp()", m_filename + ":SNP",
						genfile::string_utils::to_string( snp ),
						"entry with id " + genfile::string_utils::to_string( snp_row_id )
					) ;
				}
			}
		}
		// If we get here, the right SNP is present and we have its snp_row_id.

		for(
			std::map< std::string, std::string >::const_iterator field_i = fields.begin();
			field_i != fields.end();
		 	++field_i
 		) {
			std::string const field = field_i->first ;
			std::string const type = field_i->second ;
			get_or_create_entity( field, "" ) ;
			get_or_create_entity( type, "" ) ;
			set_relationship( field, "has_type", type ) ;
			set_relationship( field, "stored_as", "zlib_compressed_serialized" ) ;
			set_relationship( field, "is_a", "per_variant_per_sample_data" ) ;

			// Make sure we've got these fields in Entity
			statement = m_connection->get_statement(
				"SELECT id FROM Entity WHERE name == ?1"
			) ;
			statement->bind( 1, field ) ;
			statement->step() ;
			if( statement->empty() ) {
				statement = m_connection->get_statement(
					"INSERT INTO Entity( name ) VALUES( ?1 )"
				) ;
				statement->bind( 1, field ) ;
				statement->step() ;

				statement = m_connection->get_statement(
					"SELECT id FROM Entity WHERE name == ?1"
				) ;
				statement->bind( 1, field ) ;
				statement->step() ;
			}

			db::Connection::RowId field_id = statement->get_column< int >( 0 ) ;
			statement->step() ;
			if( !statement->empty() ) {
				throw genfile::DuplicateKeyError( m_filename + ":Entity", "name=\"" + field + "\"" ) ;
			}
			
			// Look to see if the data is there already.
			statement = m_connection->get_statement(
				"SELECT snp_id FROM Data WHERE snp_id == ?1 AND field_id == ?2"
			) ;
			statement->bind( 1, snp_row_id ) ;
			statement->bind( 2, field_id ) ;
			statement->step() ;
			if( statement->empty() ) {
				// std::cerr << "SNP " << snp_row_id << ", field_id " << field_id << ", statement is empty.\n" ;
				// No data already.
				// Compress the data and store it.
				std::vector< std::vector< genfile::VariantEntry > >& data = m_data ;
				data.resize( m_number_of_samples ) ;
				data_reader.get( field, data ) ;
				assert( data.size() == m_number_of_samples ) ;

				// count the data
				std::size_t count = 0 ;
				for( std::size_t i = 0; i < data.size(); ++i ) {
					count += data[i].size() ;
				}

				std::vector< char >& buffer = m_buffer ;
				m_buffer.resize( 100 + (count * 8) ) ;
				{
					char* begin = &buffer[0] ;
					char* end = begin + buffer.size() ;
					begin = genfile::write_small_integer( begin, end, uint64_t( m_number_of_samples ) ) ;
					for( std::size_t i = 0; i < data.size(); ++i ) {
						begin = genfile::write_small_integer( begin, end, data[i].size() ) ;
						for( std::size_t j = 0; j < data[i].size(); ++j ) {
							if( end < begin + 100 ) {
								std::size_t index = begin - &buffer[0] ;
								buffer.resize( buffer.size() + 1000 ) ;
								begin = &buffer[0] + index ;
								end = &buffer[0] + buffer.size() ;
							}
							begin = data[i][j].serialize( begin, end ) ;
						}
					}
					buffer.resize( begin - &buffer[0] ) ;
				}
			
				std::vector< char >& compressed_buffer = m_compressed_buffer ;
				genfile::zlib_compress( buffer, &compressed_buffer ) ;
				db::Connection::StatementPtr statement = m_connection->get_statement(
					"INSERT OR REPLACE INTO Data VALUES( ?1, ?2, ?3, ?4, ?5, ?6 )"
				) ;
				statement->bind( 1, snp_row_id ) ;
				statement->bind( 2, field_id ) ;
				statement->bind( 3, m_entities[ "unnamed cohort" ] ) ;
				statement->bind( 4, 1 ) ; // zlib-compressed, serialised.
				statement->bind( 5, uint64_t( buffer.size() ) ) ;
				statement->bind( 6, &compressed_buffer[0], compressed_buffer.size() ) ;
				statement->step() ;
			}
		}
		transaction = m_connection->get_statement(
			"COMMIT"
		) ;
		transaction->step() ;
	}
	catch( db::SQLite3Error const& e ) {
		std::cerr << "(" << e.what() << "): " << e.description() << ".\n" ;
		throw genfile::OperationFailedError( "IntensityWriter::processed_snp()", "file://" + m_filename, "SQL insert" ) ;
	}
	++m_number_of_snps_written ;
}

void IntensityWriter::end_processing_snps() {
	
}



