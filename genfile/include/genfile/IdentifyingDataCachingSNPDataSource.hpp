
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef IDENTIFYINGDATACACHINGSNPDATASOURCE_HPP
#define IDENTIFYINGDATACACHINGSNPDATASOURCE_HPP

#include "genfile/SNPDataSource.hpp"
#include "genfile/VariantIdentifyingData.hpp"
namespace genfile {
	class IdentifyingDataCachingSNPDataSource: public SNPDataSource
	{
	public:
		virtual void read_snp_identifying_data_impl( VariantIdentifyingData* variant ) = 0 ;
	private:
		void get_snp_identifying_data_impl(  VariantIdentifyingData* variant ) ;
	private:
		VariantIdentifyingData m_variant ;
	} ;
}

#endif
