#include <iostream>
#include "genfile/FileUtils.hpp"
#include "IntensityWriter.hpp"
#include "db/SQLite3Connection.hpp"
#include "db/SQLStatement.hpp"
#include "genfile/Error.hpp"
#include "genfile/zlib.hpp"

IntensityWriter::IntensityWriter( std::string const& filename ):
	m_filename( filename ),
	m_connection( new db::SQLite3Connection( filename ))
{
	setup( *m_connection ) ;
}

void IntensityWriter::setup( db::Connection& connection ) {
	db::Connection::StatementPtr statement ;
	statement = connection.get_statement(
		"CREATE TABLE IF NOT EXISTS Storage ( id INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE )"
	) ;
	statement->step() ;
	statement = connection.get_statement(
		"INSERT OR REPLACE INTO Storage VALUES( ?1, ?2 ) ;"
	) ;
	statement->bind( 1, 1 ) ;
	statement->bind( 2, "uncompressed data" ) ;
	statement->step() ;
	statement->reset() ;
	statement->bind( 1, 2 ) ;
	statement->bind( 2, "Gzip compressed data" ) ;
	statement->step() ;

	statement = connection.get_statement(
		"CREATE TABLE IF NOT EXISTS Field ( id INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE )"
	) ;
	statement->step() ;
	statement = connection.get_statement(
		"INSERT OR REPLACE INTO Field VALUES( ?1, ?2 ) ;"
	) ;
	statement->bind( 1, 0 ) ;
	statement->bind( 2, "genotypes" ) ;
	statement->step() ;
	statement->reset() ;
	statement->bind( 1, 1 ) ;
	statement->bind( 2, "intensities" ) ;
	statement->step() ;

	statement = connection.get_statement(
		"CREATE TABLE IF NOT EXISTS SNP ( id INTEGER PRIMARY KEY, rsid TEXT NOT NULL, chromosome TEXT NOT NULL, position INTEGER NOT NULL, alleleA TEXT NOT NULL, alleleB TEXT NOT NULL )"
	) ;
	statement->step() ;
	statement = connection.get_statement(
		"DELETE FROM SNP"
	) ;
	statement->step() ;
	statement = connection.get_statement(
		"CREATE TABLE IF NOT EXISTS Data ( "
			"snp_id INTEGER NOT NULL, "
			"field_id INTEGER NOT NULL, "
			"storage_id INTEGER NOT NULL, "
			"uncompressed_size INTEGER NOT NULL, "
			"data BLOB, "
			"FOREIGN KEY (field_id) REFERENCES Field( id ), "
			"FOREIGN KEY (snp_id) REFERENCES SNP( id ), "
			"FOREIGN KEY (storage_id) REFERENCES Storage( id ) "
			") ;"
	) ;
	statement->step() ;
	statement = connection.get_statement(
		"DELETE FROM Data"
	) ;
	statement->step() ;
}

void IntensityWriter::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	m_number_of_samples = number_of_samples ;
	m_number_of_snps = number_of_snps ;
	m_number_of_snps_written = 0 ;
}

void IntensityWriter::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	std::vector< std::string > fields ;
	fields.push_back( "genotypes" ) ;
	fields.push_back( "intensities" ) ;
	try {
		db::Connection::StatementPtr transaction = m_connection->get_statement(
			"BEGIN"
		) ;
		transaction->step() ;
		
		{
			db::Connection::StatementPtr statement  = m_connection->get_statement(
				"INSERT INTO SNP VALUES( ?1, ?2, ?3, ?4, ?5, ?6 )"
			) ;
			statement->bind( 1, m_number_of_snps_written ) ;
			statement->bind( 2, snp.get_rsid() ) ;
			statement->bind( 3, genfile::string_utils::to_string( snp.get_position().chromosome() )) ;
			statement->bind( 4, snp.get_position().position() ) ;
			statement->bind( 5, std::string( 1, snp.get_first_allele() ) ) ;
			statement->bind( 6, std::string( 1, snp.get_second_allele() ) ) ;
			statement->step() ;
		}
		
		for( std::size_t field_i = 0; field_i < fields.size(); ++field_i ) {
			std::string const field = fields[ field_i ] ;

			std::vector< std::vector< genfile::VariantEntry > > data( m_number_of_samples ) ;
			data_reader.get( field, genfile::VariantDataReader::set( data ) ) ;
			assert( data.size() == m_number_of_samples ) ;

			// count the data.
			std::vector< double > values ;
			values.reserve( m_number_of_samples * 3 ) ;
			for( std::size_t i = 0; i < m_number_of_samples; ++i ) {
				for( std::size_t j = 0; j < data[i].size(); ++j ) {
					if( data[i][j].is_missing() ) {
						values.push_back( std::numeric_limits< double >::quiet_NaN() ) ;
					}
					else {
						values.push_back( data[i][j].as< double >() ) ;
					}
				}
			}

			std::vector< char > buffer ;
			genfile::zlib_compress( values, &buffer ) ;
			db::Connection::StatementPtr statement = m_connection->get_statement(
				"INSERT INTO Data VALUES( ?1, ?2, ?3, ?4, ?5 )"
			) ;
			statement->bind( 1, m_number_of_snps_written ) ;
			statement->bind( 2, field_i ) ;
			statement->bind( 3, 2 ) ; // gzip compressed data.
			statement->bind( 4, buffer.size() ) ;
			statement->bind( 5, &buffer[0], buffer.size() ) ;
			statement->step() ;
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



