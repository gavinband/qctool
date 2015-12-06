
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
VariantIdentifyingData* variant	) {
		SNPIdentifyingData data ;
		m_source->get_snp_identifying_data(
			set_number_of_samples,
			set_value( data.SNPID() ),
			set_value( data.rsid() ),
			set_value( data.position().chromosome() ),
			set_value( data.position().position() ),
			set_value( data.first_allele() ),
			set_value( data.second_allele() )
		) ;
		
		Dictionary::const_iterator where = m_dictionary.find( data ) ;
		if( where != m_dictionary.end() ) {
			set_SNPID( where->second.get_SNPID() ) ;
			set_RSID( where->second.get_rsid() ) ;
			set_chromosome( where->second.get_position().chromosome() ) ;
			set_SNP_position( where->second.get_position().position() ) ;
			set_allele1( where->second.get_first_allele() ) ;
			set_allele2( where->second.get_second_allele() ) ;
		}
		else {
			set_SNPID( data.get_SNPID() ) ;
			set_RSID( data.get_rsid() ) ;
			set_chromosome( data.get_position().chromosome() ) ;
			set_SNP_position( data.get_position().position() ) ;
			set_allele1( data.get_first_allele() ) ;
			set_allele2( data.get_second_allele() ) ;
		}
	}

	VariantDataReader::UniquePtr SNPTranslatingSNPDataSource::read_variant_data_impl() {
		return m_source->read_variant_data() ;
	}

	void SNPTranslatingSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}
}
