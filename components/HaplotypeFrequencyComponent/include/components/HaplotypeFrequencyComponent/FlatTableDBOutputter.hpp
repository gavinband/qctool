
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_HAPLOTYPE_FREQUENCY_COMPONENT_FLAT_TABLE_DB_OUTPUTTER_HPP
#define QCTOOL_HAPLOTYPE_FREQUENCY_COMPONENT_FLAT_TABLE_DB_OUTPUTTER_HPP

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/bimap.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/VariantEntry.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "qcdb/Storage.hpp"
#include "qcdb/DBOutputter.hpp"

namespace haplotype_frequency_component {
	struct FlatTableDBOutputter {
		typedef std::auto_ptr< FlatTableDBOutputter > UniquePtr ;
		typedef boost::shared_ptr< FlatTableDBOutputter > SharedPtr ;
		typedef qcdb::DBOutputter::Metadata Metadata ;
		typedef qcdb::Storage::AnalysisId AnalysisId ;
		
		static UniquePtr create(
			std::string const& filename, std::string const& analysis_name, std::string const& analysis_description, Metadata const& metadata,
			std::string const& snp_match_fields = "position,alleles"
		 ) ;

		FlatTableDBOutputter(
			std::string const& filename,
			std::string const& analysis_name,
			std::string const& analysis_description,
			Metadata const& metadata,
			std::string const& snp_match_fields = "position,alleles"
		) ;

		~FlatTableDBOutputter() ;

		void set_table_name( std::string const& table_name ) ;
		
		void add_variable(
			std::string const& 
		) ;
		
		void create_new_variant_pair( genfile::VariantIdentifyingData const&, genfile::VariantIdentifyingData const& ) ;
		void store_per_variant_pair_data(
			genfile::VariantIdentifyingData const& variant1,
			genfile::VariantIdentifyingData const& variant2,
			std::string const& variable,
			genfile::VariantEntry const& value
		) ;
		
		void finalise( long options = qcdb::eCreateIndices ) ;

		AnalysisId analysis_id() const ;

	private:
		qcdb::DBOutputter m_outputter ;
		std::string m_table_name ;
		std::size_t const m_max_variants_per_block ;
		db::Connection::StatementPtr m_insert_data_sql ;
		std::vector< std::pair< genfile::VariantIdentifyingData, genfile::VariantIdentifyingData > > m_variants ;
		typedef boost::bimap< std::string, std::size_t > VariableMap ;
		VariableMap m_variables ;
		typedef std::map< std::pair< std::size_t, std::size_t >, genfile::VariantEntry > ValueMap ;
		ValueMap m_values ;

	private:
		std::string get_table_name() const ;
		void store_block() ;
		void create_schema() ;
		void store_data_for_variants(
			std::size_t const variant_i,
			db::Connection::RowId const analysis_id,
			db::Connection::RowId const variant1_id,
			db::Connection::RowId const variant2_id
		) ;
		// void create_variables() ;
	} ;
}

#endif
