
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_QCDB_FLAT_FILE_MULTIVARIANT_OUTPUTTER_HPP
#define QCTOOL_QCDB_FLAT_FILE_MULTIVARIANT_OUTPUTTER_HPP

#include <string>
#include <memory>
#include <map>
#include <vector>
#include <utility>
#include <boost/bimap.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "qcdb/StorageOptions.hpp"
#include "qcdb/MultiVariantStorage.hpp"

namespace qcdb {
	struct FlatFileMultiVariantOutputter: public MultiVariantStorage {
	public:
		typedef std::map< std::string, std::pair< std::vector< std::string >, std::string > > Metadata ;
		typedef std::auto_ptr< FlatFileMultiVariantOutputter > UniquePtr ;
		static UniquePtr create(
			std::string const& filename,
			std::size_t const number_of_key_entries,
			std::string const& analysis_name,
			Metadata const& metadata
		) ;

	public:
		FlatFileMultiVariantOutputter(
			std::string const& filename,
			std::size_t const number_of_key_entries,
			std::string const& analysis_name,
			Metadata const& metadata
		) ;
		~FlatFileMultiVariantOutputter() ;

		void set_variant_names( std::vector< std::string > const& names ) ;
		void add_variable( std::string const& ) ;

		void create_new_key( Key const& key ) ;
		void store_data_for_key(
			Key const& key,
			std::string const& value_name,
			genfile::VariantEntry const& value
		) ;

		void finalise( long options = eCreateIndices ) ;

		AnalysisId analysis_id() const ;

	private:
		std::string const m_filename ;
		std::vector< std::string > m_key_entry_names ;
		std::string const m_analysis_name ;
		Metadata const m_metadata ;
		std::size_t const m_max_keys_per_block ;
		statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
		std::vector< Key > m_keys ;
		typedef boost::bimap< std::string, std::size_t > VariableMap ;
		VariableMap m_variables ;
		typedef std::map< std::pair< std::size_t, std::size_t >, genfile::VariantEntry > ValueMap ;
		ValueMap m_values ;
	
	private:
		void store_block() ;
		std::string format_metadata() const ;
		void create_columns() ;
	} ;
}

#endif
