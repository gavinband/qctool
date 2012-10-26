
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_STORAGE_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_STORAGE_HPP

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include "genfile/SNPIdentifyingData2.hpp"

namespace snp_summary_component {
	struct Storage {
		typedef std::auto_ptr< Storage > UniquePtr ;
		typedef boost::shared_ptr< Storage > SharedPtr ;
		
		virtual ~Storage() {} ;

		virtual void store_per_variant_data(
			genfile::SNPIdentifyingData2 const& snp,
			std::string const& variable,
			genfile::VariantEntry const& value
		) = 0 ;
		
		virtual void finalise() {} ;
	} ;
}

#endif

