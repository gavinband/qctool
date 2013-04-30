
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_QCDB_FLAT_FILE_OUTPUTTER_HPP
#define QCTOOL_QCDB_FLAT_FILE_OUTPUTTER_HPP

#include <string>
#include <memory>
#include <map>
#include <vector>
#include <utility>
#include <boost/bimap.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/SNPIdentifyingData2.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "qcdb/Storage.hpp"
#include "qcdb/StorageOptions.hpp"

namespace qcdb {
	struct FlatFileOutputter: public Storage {
		typedef std::map< std::string, std::pair< std::vector< std::string >, std::string > > Metadata ;
		static UniquePtr create( std::string const& filename, std::string const& analysis_name, Metadata const& metadata ) ;
		static SharedPtr create_shared( std::string const& filename, std::string const& analysis_name, Metadata const& metadata ) ;

		FlatFileOutputter( std::string const& filename, std::string const& analysis_name, Metadata const& metadata ) ;
		~FlatFileOutputter() ;

		void add_variable( std::string const& ) ;

		void create_new_variant( genfile::SNPIdentifyingData2 const& ) ;
		void store_per_variant_data(
			genfile::SNPIdentifyingData2 const& snp,
			std::string const& value_name,
			genfile::VariantEntry const& value
		) ;

		void finalise( long options = eCreateIndices ) ;

	private:
		std::string const m_filename ;
		std::string const m_analysis_name ;
		Metadata const m_metadata ;
		std::size_t const m_max_snps_per_block ;
		statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
		std::vector< genfile::SNPIdentifyingData2 > m_snps ;
		typedef boost::bimap< std::string, std::size_t > VariableMap ;
		VariableMap m_variables ;
		typedef std::map< std::pair< std::size_t, std::size_t >, genfile::VariantEntry > ValueMap ;
		ValueMap m_values ;
	
	private:
		void store_block() ;
		std::string format_metadata() const ;
	} ;
}

#endif
