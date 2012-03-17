#ifndef QCTOOL_SAMPLE_SUMMARY_COMPONENT_DB_OUTPUTTER_HPP
#define QCTOOL_SAMPLE_SUMMARY_COMPONENT_DB_OUTPUTTER_HPP

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

namespace sample_stats {
	struct DBOutputter: public qcdb::DBOutputter {
		typedef std::auto_ptr< DBOutputter > UniquePtr ;
		typedef boost::shared_ptr< DBOutputter > SharedPtr ;

		static UniquePtr create( std::string const& filename, std::string const& cohort_name, Metadata const& metadata, genfile::CohortIndividualSource const& samples ) ;
		static SharedPtr create_shared( std::string const& filename, std::string const& cohort_name, Metadata const& metadata, genfile::CohortIndividualSource const& samples ) ;

		DBOutputter( std::string const& filename, std::string const& cohort_name, Metadata const& metadata, genfile::CohortIndividualSource const& samples ) ;
		~DBOutputter() ;

		void operator()(
			std::string const& computation_name,
			std::size_t sample,
			std::string const& variable,
			std::string const& description,
			genfile::VariantEntry const& value
		) ;

	private:
		genfile::CohortIndividualSource const& m_samples ;
		std::size_t const m_max_transaction_count ;
		db::Connection::StatementPtr m_find_sample_statement ;
		db::Connection::StatementPtr m_insert_sample_statement ;
		db::Connection::StatementPtr m_find_sampledata_statement ;
		db::Connection::StatementPtr m_insert_sampledata_statement ;

		typedef std::vector< boost::tuple< std::string, std::size_t, std::string, std::string, genfile::VariantEntry > > Data ;
		Data m_data ;

	private:
		void construct_statements() ;
		void store_samples( genfile::CohortIndividualSource const& samples ) ;
		void store_sample( genfile::CohortIndividualSource const& samples, std::size_t sample ) ;
		db::Connection::RowId get_or_create_sample( genfile::VariantEntry const& identifier ) const ;
		void write_data( Data const& data ) ;
		void write_data(
			std::size_t sample,
			std::string const& variable,
			std::string const& description,
			genfile::VariantEntry const& value
		) ;
		
	} ;
}

#endif
