#include <iostream>
#include <boost/bind.hpp>
#include "genfile/FileUtils.hpp"
#include "db/SQLite3Connection.hpp"
#include "db/SQLStatement.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"
#include "genfile/zlib.hpp"
#include "genfile/endianness_utils.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "DataStore.hpp"
#include "VCDBWriter.hpp"

void VCDBWriter::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Database options" ) ;
	options[ "-write-db" ]
		.set_description( "Write data from the input files to the given database." )
		.set_takes_single_value() ;
	options[ "-write-db-fields" ]
		.set_description( "Specify fields from the input files to write." )
		.set_takes_single_value()
		.set_default_value( "*" ) ;
	options.option_implies_option( "-write-db-fields", "-write-db" ) ;
}

VCDBWriter::UniquePtr VCDBWriter::create( appcontext::OptionProcessor const& options ) {
	std::string const filename = options.get< std::string >( "-write-db" ) ;
	std::vector< std::string > const fields = genfile::string_utils::split_and_strip(
		options.get_value< std::string >( "-write-db-fields" ), ",", " \t"
	) ;
	VCDBWriter::UniquePtr result(
		new VCDBWriter(
			DataStore::create(
				"sqlite3://" + filename
			),
			fields
		)
	) ;
	return result ;
}

VCDBWriter::VCDBWriter( DataStore::UniquePtr store, std::vector< std::string > const& fields ):
	m_store( store ),
	m_fields( fields.begin(), fields.end() ),
	m_write_all_fields( m_fields.size() == 1 && (m_fields.begin()->compare( "*" ) == 0 ) )
{
	setup() ;
}

void VCDBWriter::setup() {
	m_storage_id = m_store->get_or_create_entity( "qctool_serialized", "serialized using qctool." ) ;
	m_store->set_relationship( "qctool_serialized", "is_a", "storage_type" ) ;
	m_store->get_or_create_entity( "per_variant_per_sample_data", "Indicates data pertains to each sample at a variant" ) ;	
	m_cohort_id = m_store->get_or_create_entity( "unnamed cohort", "An unnamed cohort" ) ;	
}

db::Connection::RowId VCDBWriter::get_or_create_field( std::string const& field, std::string const& type ) {
	EntityCache::const_iterator where = m_entity_cache.find( field ) ;
	if( where != m_entity_cache.end() ) {
		return where->second ;
	}
	else {
		DataStore::EntityId field_id = m_store->get_or_create_entity( field, "" ) ;
		m_store->get_or_create_entity( type, "" ) ;
		m_store->set_relationship( field, "has_type", type ) ;
		m_store->set_relationship( field, "stored_as", "qctool_serialized" ) ;
		m_store->set_relationship( field, "is_a", "per_variant_per_sample_data" ) ;
		m_entity_cache[ field ] = field_id ;
		return field_id ;
	}
}

void VCDBWriter::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	m_number_of_samples = number_of_samples ;
	m_number_of_snps = number_of_snps ;
	m_number_of_snps_written = 0 ;
	m_transaction_count = std::numeric_limits< std::size_t >::max() ;
	m_transaction_limit = 100 ;
}

void VCDBWriter::set_SNPs( std::vector< genfile::SNPIdentifyingData > const& snps ) {
	std::size_t const CHUNK_SIZE = 1000 ;
	for( std::size_t chunk = 0; chunk < ( snps.size() / CHUNK_SIZE ); ++chunk ) {
		DataStore::TransactionPtr transaction = m_store->open_transaction() ;
		for( std::size_t i = chunk * CHUNK_SIZE; i < std::min( (chunk+1) * CHUNK_SIZE, snps.size() ); ++i ) {
			m_store->get_or_create_SNP( snps[i] ) ;
		}
	}
}

namespace impl {
	void add_spec_to_map( std::map< std::string, std::string >* map, std::string const& name, std::string const& type ) {
		(*map)[ name ] = type ;
	}
}

void VCDBWriter::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	std::map< std::string, std::string > fields ;
	data_reader.get_supported_specs( boost::bind( &impl::add_spec_to_map, &fields, _1, _2 )) ;
	
	if( m_transaction_count >= m_transaction_limit ) {
		// Commit the current transaction and start a new one.
		DataStore::TransactionPtr transaction = m_store->open_transaction() ;
		m_transaction_count = 0 ;
	}
	DataStore::EntityId snp_id  = m_store->get_or_create_SNP( snp ) ;

	for(
		std::map< std::string, std::string >::const_iterator field_i = fields.begin();
		field_i != fields.end();
	 	++field_i
		) {
		if( m_write_all_fields || m_fields.find( field_i->first ) != m_fields.end() ) {
			store_field( snp_id, field_i->first, field_i->second, data_reader ) ;
			++m_transaction_count ;
		}
	}
}

void VCDBWriter::store_field( DataStore::EntityId const snp_id, std::string const& field, std::string const& type, genfile::VariantDataReader& data_reader ) {
	try {
		unsafe_store_field( snp_id, field, type, data_reader ) ;
	}
	catch( db::StatementPreparationError const& e ) {
		std::cerr << "(" << e.what() << "): " << e.description() << ": preparing statement \"" + e.SQL() + "\".\n" ;
		throw genfile::OperationFailedError( "VCDBWriter::processed_snp()", m_store->get_spec(), "SQL insert" ) ;
	}
	catch( db::Error const& e ) {
		std::cerr << "(" << e.what() << "): " << e.description() << ".\n" ;
		throw genfile::OperationFailedError( "VCDBWriter::processed_snp()", m_store->get_spec(), "SQL insert" ) ;
	}
}
	
void VCDBWriter::unsafe_store_field( DataStore::EntityId const snp_id, std::string const& field, std::string const& type, genfile::VariantDataReader& data_reader ) {
	DataStore::EntityId const field_id = get_or_create_field( field, type ) ;

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

	m_store->store_per_variant_data(
		snp_id,
		field_id,
		m_cohort_id,
		m_storage_id,
		&buffer[0],
		&buffer[0] + buffer.size()
	) ;
	++m_number_of_snps_written ;
}

void VCDBWriter::end_processing_snps() {
	
}



