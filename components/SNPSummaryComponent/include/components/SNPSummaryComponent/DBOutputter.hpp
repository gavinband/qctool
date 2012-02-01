#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_DB_OUTPUTTER_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_DB_OUTPUTTER_HPP

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"

namespace impl {
	struct DBOutputter {
		typedef std::auto_ptr< DBOutputter > UniquePtr ;
		typedef boost::shared_ptr< DBOutputter > SharedPtr ;

		static UniquePtr create( std::string const& filename, std::string const& cohort_name, std::string const& data_source_spec, std::string const& exclusions_name ) ;
		static SharedPtr create_shared( std::string const& filename, std::string const& cohort_name, std::string const& data_source_spec, std::string const& exclusions_name ) ;

		DBOutputter( std::string const& filename, std::string const& cohort_name, std::string const& data_source_spec, std::string const& exclusions_name ) ;
		~DBOutputter() ;

		void operator()(
			std::size_t index,
			genfile::SNPIdentifyingData const& snp,
			std::string const& computation_name,
			std::string const& variable,
			genfile::VariantEntry const& value
		) ;

	private:
		db::Connection::UniquePtr m_connection ;
		std::size_t const m_max_transaction_count ;
		std::string const m_cohort_name ;
		std::string const m_source_spec ;
		std::string const m_exclusions_name ;

		db::Connection::StatementPtr m_find_variant_statement ;
		db::Connection::StatementPtr m_insert_variant_statement ;
		db::Connection::StatementPtr m_find_entity_statement ;
		db::Connection::StatementPtr m_find_entity_data_statement ;
		db::Connection::StatementPtr m_insert_entity_statement ;
		db::Connection::StatementPtr m_insert_entity_data_statement ;
		db::Connection::StatementPtr m_insert_summarydata_statement ;
		db::Connection::RowId m_analysis_id ;
		typedef std::vector< boost::tuple< genfile::SNPIdentifyingData, std::string, genfile::VariantEntry > > Data ;
		Data m_data ;

	private:
		void construct_statements() ;
		void reset_statements() ;
		void write_data( Data const& data ) ;

		db::Connection::RowId get_or_create_snp( genfile::SNPIdentifyingData const& snp ) const ;
		db::Connection::RowId get_or_create_variable( std::string const& name ) const ;
		db::Connection::RowId get_or_create_entity( std::string const& name ) const ;
		db::Connection::RowId get_or_create_entity( std::string const& name, std::string const& description ) const ;
		db::Connection::RowId get_or_create_entity_data( db::Connection::RowId const entity_id, db::Connection::RowId const variable_id, genfile::VariantEntry const& value ) const ;

		void store_data(
			genfile::SNPIdentifyingData const& snp,
			std::string const& variable,
			genfile::VariantEntry const& value
		) ;
	} ;
}

#endif
