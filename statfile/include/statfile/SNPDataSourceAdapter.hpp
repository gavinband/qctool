
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef STATFILE_SNPDATASOURCE_ADAPTER_HPP
#define STATFILE_SNPDATASOURCE_ADAPTER_HPP

#include "genfile/SNPDataSource.hpp"
#include "genfile/Error.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"

namespace statfile {
	struct SNPDataSourceAdapter: public genfile::SNPDataSource {
		SNPDataSourceAdapter( BuiltInTypeStatSource::UniquePtr source ) ;

		operator bool() const { return *m_source ; }
		unsigned int number_of_samples() const { return 0 ; }
		OptionalSnpCount total_number_of_snps() const { return m_source->number_of_rows() ; }

		std::string get_source_spec() const ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
		void get_snp_identifying_data_impl( 
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			ChromosomeSetter const& set_chromosome,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) ;	

		void read_snp_probability_data_impl(
			GenotypeProbabilitySetter const& ,
			std::string const& genotype_field
		) ;

		genfile::VariantDataReader::UniquePtr read_variant_data_impl() ;
		void ignore_snp_probability_data_impl() ;
		void reset_to_start_impl() ;
		
	private:
		BuiltInTypeStatSource::UniquePtr m_source ;
		genfile::SNPIdentifyingData m_snp ;
	private:
		void setup() ;
	} ;

}

#endif
