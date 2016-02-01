
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SNPTRANSLATING_SNP_DATA_SOURCE_HPP
#define GENFILE_SNPTRANSLATING_SNP_DATA_SOURCE_HPP

#include <map>
#include <string>
#include <memory>

#include "genfile/SNPDataSource.hpp"
#include "genfile/VariantIdentifyingData.hpp"

namespace genfile {

	class SNPTranslatingSNPDataSource: public SNPDataSource
	{
	public:
		typedef std::auto_ptr< SNPTranslatingSNPDataSource > UniquePtr ;
		typedef std::map< VariantIdentifyingData, VariantIdentifyingData > Dictionary ;
		
		static UniquePtr create(
			SNPDataSource::UniquePtr source,
			Dictionary const& dict
		) ;

	public:
		SNPTranslatingSNPDataSource(
			SNPDataSource::UniquePtr source,
			Dictionary const& dict
		) ;
	
		operator bool() const {
			return  (*m_source) ;
		}

		Metadata get_metadata() const {
			return m_source->get_metadata() ;
		}
		
		unsigned int number_of_samples() const {
			return m_source->number_of_samples() ;
		}

		OptionalSnpCount total_number_of_snps() const {
			return m_source->total_number_of_snps() ;
		}

		std::string get_source_spec() const {
			return "snp-translated:" + m_source->get_source_spec() ;
		}

		SNPDataSource const& get_parent_source() const {
			return *m_source ;
		}

		SNPDataSource const& get_base_source() const {
			return m_source->get_base_source() ;
		}

	private:
		SNPDataSource::UniquePtr m_source ;
		Dictionary m_dictionary ;

	private:
		void reset_to_start_impl() {
			m_source->reset_to_start() ;
		}

		void get_snp_identifying_data_impl( 
			VariantIdentifyingData* variant
		) ;

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() ;
	} ;
}

#endif
