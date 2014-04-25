
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SAMPLE_SUMMARY_COMPONENT_SAMPLE_STORAGE_HPP
#define QCTOOL_SAMPLE_SUMMARY_COMPONENT_SAMPLE_STORAGE_HPP

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include "genfile/VariantEntry.hpp"
#include "qcdb/StorageOptions.hpp"

namespace sample_stats {
	struct SampleStorage {
		typedef std::auto_ptr< SampleStorage > UniquePtr ;
		typedef boost::shared_ptr< SampleStorage > SharedPtr ;

		virtual ~SampleStorage() {} ;
		virtual void operator()(
			std::string const& computation_name,
			std::size_t sample,
			std::string const& variable,
			std::string const& description,
			genfile::VariantEntry const& value
		) = 0 ;
		
		virtual void finalise( long options = qcdb::eCreateIndices ) = 0 ;
	} ;
}

#endif
