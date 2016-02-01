
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/IdentifyingDataCachingSNPDataSource.hpp"

namespace genfile {
	void IdentifyingDataCachingSNPDataSource::get_snp_identifying_data_impl( VariantIdentifyingData* variant ) {
		if( state() != e_HaveReadIdentifyingData ) {
			read_snp_identifying_data_impl( &m_variant ) ;
		}
		*variant = m_variant ;
	}
}
