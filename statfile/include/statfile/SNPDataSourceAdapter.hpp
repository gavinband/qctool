
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef STATFILE_SNPDATASOURCE_ADAPTER_HPP
#define STATFILE_SNPDATASOURCE_ADAPTER_HPP

#include "genfile/SNPDataSource.hpp"
#include "genfile/Error.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"

namespace statfile {
	struct SNPDataSourceAdapter: public genfile::SNPDataSource {
		SNPDataSourceAdapter( BuiltInTypeStatSource::UniquePtr source ) ;

		operator bool() const { return *m_source ; }
		Metadata get_metadata() const { return Metadata() ; }
		unsigned int number_of_samples() const { return 0 ; }
		OptionalSnpCount total_number_of_snps() const { return m_source->number_of_rows() ; }

		std::string get_source_spec() const ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;

		void get_snp_identifying_data_impl( genfile::VariantIdentifyingData* result ) ;

		void read_snp_probability_data_impl(
			GenotypeProbabilitySetter const& ,
			std::string const& genotype_field
		) ;

		genfile::VariantDataReader::UniquePtr read_variant_data_impl() ;
		void ignore_snp_probability_data_impl() ;
		void reset_to_start_impl() ;
		
	private:
		BuiltInTypeStatSource::UniquePtr m_source ;
		genfile::VariantIdentifyingData m_snp ;
	private:
		void setup() ;
	} ;

}

#endif
