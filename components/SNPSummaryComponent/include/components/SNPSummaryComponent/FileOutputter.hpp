
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_FILE_OUTPUTTER_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_FILE_OUTPUTTER_HPP

#include <string>
#include <memory>
#include "genfile/VariantEntry.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "components/SNPSummaryComponent/Storage.hpp"

namespace snp_summary_component {
	struct FileOutputter: public Storage {
		static UniquePtr create( std::string const& filename ) { return UniquePtr( new FileOutputter( filename ) ) ; }
		static SharedPtr create_shared( std::string const& filename ) { return SharedPtr( new FileOutputter( filename ) ) ; }

		FileOutputter( std::string const& filename ):
			m_filename( filename ),
			m_sink( statfile::BuiltInTypeStatSink::open( filename ))
		{
			(*m_sink) | "SNPID" | "rsid" | "chromosome" | "position" | "alleleA" | "alleleB" | "variable" | "value" ;
		}

		void store_per_variant_data(
			genfile::VariantIdentifyingData const& snp,
			std::string const& value_name,
			genfile::VariantEntry const& value
		) {
			std::vector< genfile::string_utils::slice > const alternate_ids = snp.get_identifiers() ;
			std::string SNPID ;
			for( std::size_t i = 0; i < alternate_ids.size(); ++i ) {
				SNPID += ( i > 0 ? "," : "" ) + std::string( alternate_ids[i] ) ;
			}
			(*m_sink)
				<< SNPID
				<< snp.get_rsid()
				<< std::string( snp.get_position().chromosome() )
				<< snp.get_position().position()
				<< snp.get_first_allele()
				<< snp.get_second_allele()
				<< value_name ;
			if( value == value ) {
				(*m_sink) << value ;
			}
			else {
				(*m_sink) << "NA" ;
			}
			(*m_sink) << statfile::end_row() ;
			;
		}

	private:
		std::string const m_filename ;
		statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
	} ;
}

#endif
