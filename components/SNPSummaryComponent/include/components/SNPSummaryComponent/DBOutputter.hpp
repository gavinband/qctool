#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_DB_OUTPUTTER_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_DB_OUTPUTTER_HPP

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/VariantEntry.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "qcdb/DBOutputter.hpp"

namespace impl {
	struct DBOutputter: public qcdb::DBOutputter {
		typedef std::auto_ptr< DBOutputter > UniquePtr ;
		typedef boost::shared_ptr< DBOutputter > SharedPtr ;
		
		static UniquePtr create( std::string const& filename, std::string const& cohort_name, Metadata const& metadata ) ;
		static SharedPtr create_shared( std::string const& filename, std::string const& cohort_name, Metadata const& metadata ) ;

		DBOutputter(
			std::string const& filename,
			std::string const& cohort_name,
			Metadata const& metadata
		) ;
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
		db::Connection::StatementPtr m_find_variant_statement ;
		db::Connection::StatementPtr m_insert_variant_statement ;
		typedef std::vector< boost::tuple< genfile::SNPIdentifyingData, std::string, genfile::VariantEntry > > Data ;
		Data m_data ;

	private:
		void construct_statements() ;
		void write_data( Data const& data ) ;

		db::Connection::RowId get_or_create_snp( genfile::SNPIdentifyingData const& snp ) const ;

		void store_data(
			genfile::SNPIdentifyingData const& snp,
			std::string const& variable,
			genfile::VariantEntry const& value
		) ;
	} ;
}

#endif
