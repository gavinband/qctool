#ifndef VCDB_DATA_STORE_HPP
#define VCDB_DATA_STORE_HPP

#include <map>
#include <string>
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "DataStore.hpp"

struct VCDBDataStore: public DataStore 
{
	~VCDBDataStore() throw() {}
	VCDBDataStore( db::Connection::UniquePtr connection ) ;
	std::string get_spec() const ;
	int64_t get_or_create_SNP( genfile::SNPIdentifyingData const& ) ;
	int64_t get_or_create_entity( std::string name, std::string description ) ;
	int64_t get_entity( std::string const& name ) const ;
	void set_relationship( std::string const& left, std::string const& relation, std::string const& right ) const ;
	void store_per_variant_data( int64_t snp_id, std::string const& field, std::string const& cohort, std::string const& storage, char const* buffer, char const* const end ) ;
	void store_per_variant_summary_data( int64_t snp_id, std::string const& tool, std::string const& field, std::string const& storage, genfile::VariantEntry const& value ) ;
	void get_entities_by_relation( std::string const& relationship, std::string const& related_entity, boost::function< void ( db::Connection::RowId, std::string const& ) > callback ) ;
	Transaction::UniquePtr open_transaction() ;
private:
	db::Connection::UniquePtr m_connection ;
	std::vector< char > m_compression_buffer ;
	std::string m_db_version ;
	db::Connection::StatementPtr m_store_variant_data_statement ;
	db::Connection::RowId m_zlib_compression_id ;
	db::Connection::RowId m_no_compression_id ;
	std::vector< std::pair< genfile::SNPIdentifyingData, db::Connection::RowId > > m_variants ;

private:
	std::string get_db_version( db::Connection& connection ) const ;
	void prepare_new_db( db::Connection& connection, std::string const& version ) ;
	void prepare_existing_db( db::Connection& connection, std::string const& version ) ;
	void prepare_db( db::Connection& connection, std::string const& version, bool new_db ) ;
	void create_indices( db::Connection& connection ) ;
	void setup_db( db::Connection& connection ) ;
	void load_variants( db::Connection& connection ) ;
	void set_relationship( db::Connection::RowId left, db::Connection::RowId relation, db::Connection::RowId right ) const ;
	void store_per_variant_data( int64_t snp_id, int64_t field_id, int64_t cohort_id, int64_t storage_id, char const* buffer, char const* const end ) ;
	void get_entities_by_relation( int64_t relationship, int64_t related_entity, boost::function< void ( db::Connection::RowId, std::string const& ) > callback ) ;
	db::Connection::StatementPtr get_store_variant_data_statement( std::string const& db_version) ;

	struct VCDBTransaction:public DataStore::Transaction
	{
		typedef std::auto_ptr< VCDBTransaction > UniquePtr ;
		VCDBTransaction( db::Connection& connection ) ;
		~VCDBTransaction() ;
	private:
		db::Connection* m_connection ;
		bool m_rollback ;
	} ;
} ;

#endif

