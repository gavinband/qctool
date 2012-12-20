
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_DB_OUTPUTTER_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_DB_OUTPUTTER_HPP

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include "genfile/SNPIdentifyingData2.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/VariantEntry.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "qcdb/DBOutputter.hpp"
#include "qcdb/Storage.hpp"

namespace snp_summary_component {
	struct DBOutputter: public qcdb::Storage {
		typedef qcdb::DBOutputter::Metadata Metadata ;
		typedef qcdb::Storage Storage ;
		static Storage::UniquePtr create( std::string const& filename, std::string const& cohort_name, Metadata const& metadata ) ;
		static Storage::SharedPtr create_shared( std::string const& filename, std::string const& cohort_name, Metadata const& metadata ) ;

		DBOutputter(
			std::string const& filename,
			std::string const& cohort_name,
			Metadata const& metadata
		) ;

		~DBOutputter() ;

		void store_per_variant_data(
			genfile::SNPIdentifyingData2 const& snp,
			std::string const& variable,
			genfile::VariantEntry const& value
		) ;
		
		void finalise() ;

	private:
		qcdb::DBOutputter m_outputter ;
		std::size_t const m_max_transaction_count ;
		db::Connection::RowId m_variable_id ;
		typedef std::vector< boost::tuple< genfile::SNPIdentifyingData2, std::string, genfile::VariantEntry > > Data ;
		Data m_data ;

	private:
		void write_data( Data const& data ) ;

		void store_data(
			db::Connection::RowId const snp_id,
			std::string const& variable,
			genfile::VariantEntry const& value
		) ;
	} ;
}

#endif
