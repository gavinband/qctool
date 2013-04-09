#ifndef COMPONENTS_HAPLOTYPEFREQUENCYCOMPONENT_HAPLOTYPEFREQUENCYDBOUTPUTTER_HPP
#define COMPONENTS_HAPLOTYPEFREQUENCYCOMPONENT_HAPLOTYPEFREQUENCYDBOUTPUTTER_HPP

#include <utility>
#include <string>
#include <genfile/VariantEntry.hpp>
#include <boost/shared_ptr.hpp>
#include "qcdb/DBOutputter.hpp"

namespace haplotype_frequency_component {
	struct DBOutputter: public qcdb::DBOutputter {

		typedef std::auto_ptr< DBOutputter > UniquePtr ;
		typedef boost::shared_ptr< DBOutputter > SharedPtr ;

		static UniquePtr create( std::string const& filename, std::string const& analysis_name, std::string const& analysis_description, Metadata const& metadata ) ;
		static SharedPtr create_shared( std::string const& filename, std::string const& analysis_name, std::string const& analysis_description, Metadata const& metadata ) ;

		DBOutputter( std::string const& filename, std::string const& analysis_name, std::string const& analysis_description, Metadata const& metadata ) ;

		~DBOutputter() ;

		void operator()(
			std::string const& cohort,
			genfile::SNPIdentifyingData const& source_snp,
			genfile::SNPIdentifyingData const& target_snp,
			std::string const& variable,
			genfile::VariantEntry const& value
		) ;

	private:
		std::size_t const m_max_transaction_count ;
		db::Connection::RowId m_variable_class_id ;
		db::Connection::StatementPtr m_insert_summarydata_statement ;

		typedef std::vector< boost::tuple< std::string, genfile::SNPIdentifyingData, genfile::SNPIdentifyingData, std::string, genfile::VariantEntry > > Data ;
		Data m_data ;

	private:
		void construct_statements() ;
		void write_data( Data const& data ) ;
		void reset_statements() ;
		void store_comparison(
			db::Connection::RowId const variant1_id,
			db::Connection::RowId const variant2_id,
			db::Connection::RowId const analysis_id,
			db::Connection::RowId const variable_id,
			genfile::VariantEntry const& value
		) ;
	} ;
}

#endif
