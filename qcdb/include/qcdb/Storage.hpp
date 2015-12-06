
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_QCDB_STORAGE_HPP
#define QCTOOL_QCDB_STORAGE_HPP

#include <stdint.h>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "qcdb/StorageOptions.hpp"

namespace qcdb {
	struct Storage {
		typedef std::auto_ptr< Storage > UniquePtr ;
		typedef boost::shared_ptr< Storage > SharedPtr ;
		
		virtual ~Storage() {} ;

		virtual void add_variable(
			std::string const& 
		) = 0 ;

		virtual void create_new_variant( genfile::VariantIdentifyingData const& ) = 0 ;
		virtual void store_per_variant_data(
			genfile::VariantIdentifyingData const& snp,
			std::string const& variable,
			genfile::VariantEntry const& value
		) = 0 ;
		
		virtual void finalise( long options = eCreateIndices ) {} ;
		
		typedef int64_t AnalysisId ;
		virtual AnalysisId analysis_id() const = 0 ;
	} ;
}

#endif

