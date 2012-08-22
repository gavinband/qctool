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
	} ;
}

#endif

