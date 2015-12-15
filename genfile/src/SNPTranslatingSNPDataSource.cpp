
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <map>
#include <string>
#include <memory>

#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPTranslatingSNPDataSource.hpp"
#include "genfile/get_set.hpp"

namespace genfile {
	SNPTranslatingSNPDataSource::UniquePtr SNPTranslatingSNPDataSource::create(
		SNPDataSource::UniquePtr source,
		Dictionary const& dict
	) {
		return UniquePtr(
			new SNPTranslatingSNPDataSource( source, dict )
		) ;
	}

	SNPTranslatingSNPDataSource::SNPTranslatingSNPDataSource(
		SNPDataSource::UniquePtr source,
		Dictionary const& dict
	):
		m_source( source ),
		m_dictionary( dict )
	{
		assert( m_source.get() ) ;
	}
		
	void SNPTranslatingSNPDataSource::get_snp_identifying_data_impl( 
		VariantIdentifyingData* variant
	) {
		VariantIdentifyingData data ;
		m_source->get_snp_identifying_data( &data ) ;

		Dictionary::const_iterator where = m_dictionary.find( data ) ;
		if( where != m_dictionary.end() ) {
			*variant = where->second ;
		}
		else {
			*variant = data ;
		}
	}

	VariantDataReader::UniquePtr SNPTranslatingSNPDataSource::read_variant_data_impl() {
		return m_source->read_variant_data() ;
	}

	void SNPTranslatingSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}
}
