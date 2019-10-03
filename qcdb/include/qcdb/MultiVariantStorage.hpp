
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_QCDB_MULTIVARIANT_STORAGE_HPP
#define QCTOOL_QCDB_MULTIVARIANT_STORAGE_HPP

#include <string>
#include <memory>
#include <map>
#include <stdint.h>
#include <boost/shared_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "qcdb/StorageOptions.hpp"
#include "qcdb/Storage.hpp"
#include "qcdb/DBOutputter.hpp"

namespace qcdb {
	struct MultiVariantStorage {
		typedef std::auto_ptr< MultiVariantStorage > UniquePtr ;
		typedef Storage::Metadata Metadata ;
		typedef Storage::AnalysisId AnalysisId ;
		typedef std::vector< genfile::VariantIdentifyingData > Key ;
		
		static UniquePtr create(
			std::string const& filename,
			std::size_t const number_of_key_fields,
			std::string const& analysis_name,
			std::string const& analysis_chunk,
			Metadata const& metadata,
			std::string const& compare_by = "position,alleles",
			boost::optional< genfile::db::Connection::RowId > analysis_id = boost::optional< genfile::db::Connection::RowId >()
		) ;
		
		virtual ~MultiVariantStorage() {} ;

		virtual void set_variant_names( std::vector< std::string > const& names ) = 0 ;

		virtual void add_variable(
			std::string const& 
		) = 0 ;
		virtual void create_new_key(
			Key const& key
		) = 0 ;
		virtual void store_data_for_key(
			Key const& key,
			std::string const& variable,
			genfile::VariantEntry const& value
		) = 0 ;
		
		virtual void finalise( long options = eCreateIndices ) {} ;
		
		virtual AnalysisId analysis_id() const = 0 ;
	} ;
}

#endif

