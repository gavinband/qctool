#ifndef QCTOOL_SAMPLE_SUMMARY_COMPONENT_SAMPLE_DB_OUTPUTTER_HPP
#define QCTOOL_SAMPLE_SUMMARY_COMPONENT_SAMPLE_DB_OUTPUTTER_HPP

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "qcdb/DBOutputter.hpp"

namespace impl {
	struct SampleDBOutputter: public qcdb::DBOutputter {
		typedef std::auto_ptr< SampleDBOutputter > UniquePtr ;
		typedef boost::shared_ptr< SampleDBOutputter > SharedPtr ;

		static UniquePtr create( std::string const& filename, std::string const& cohort_name, std::string const& data_source_spec, std::string const& exclusions_name ) ;
		static SharedPtr create_shared( std::string const& filename, std::string const& cohort_name, std::string const& data_source_spec, std::string const& exclusions_name ) ;

		SampleDBOutputter( std::string const& filename, std::string const& cohort_name, std::string const& data_source_spec, std::string const& exclusions_name ) ;
		~SampleDBOutputter() ;

		void operator()(
			std::size_t index,
			genfile::SNPIdentifyingData const& snp,
			std::string const& computation_name,
			std::string const& variable,
			genfile::VariantEntry const& value
		) ;

		void store_samples( genfile::CohortIndividualSource const& samples ) ;

	private:
		db::Connection::UniquePtr m_connection ;
		std::size_t const m_max_transaction_count ;
		std::string const m_cohort_name ;
		std::string const m_source_spec ;
		std::string const m_exclusions_name ;

		db::Connection::StatementPtr m_find_sample_statement ;
		db::Connection::StatementPtr m_insert_sample_statement ;
		db::Connection::StatementPtr m_find_entity_statement ;
		db::Connection::StatementPtr m_find_entity_with_description_statement ;
		db::Connection::StatementPtr m_find_entity_data_statement ;
		db::Connection::StatementPtr m_insert_entity_statement ;
		db::Connection::StatementPtr m_insert_entity_data_statement ;
		db::Connection::StatementPtr m_insert_sampledata_statement ;
		db::Connection::RowId m_analysis_id ;
		typedef std::vector< boost::tuple< genfile::SNPIdentifyingData, std::string, genfile::VariantEntry > > Data ;
		Data m_data ;

	private:
		void construct_statements() ;
		void reset_statements() ;
		void store_sample( genfile::CohortIndividualSource const& samples, std::size_t sample ) ;

		db::Connection::RowId get_or_create_sample( genfile::VariantEntry const& identifier ) const ;
		db::Connection::RowId get_or_create_variable( std::string const& name, std::string const& description ) const ;
		db::Connection::RowId get_or_create_entity( std::string const& name ) const ;
		db::Connection::RowId get_or_create_entity( std::string const& name, std::string const& description ) const ;
		db::Connection::RowId get_or_create_entity_data( db::Connection::RowId const entity_id, db::Connection::RowId const variable_id, genfile::VariantEntry const& value ) const ;
	} ;
}

#endif
