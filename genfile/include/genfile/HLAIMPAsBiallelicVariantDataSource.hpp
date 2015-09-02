
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef HLAIMPSNPDATASOURCE_HPP
#define HLAIMPSNPDATASOURCE_HPP

#include <iostream>
#include <string>
#include <map>
#include <boost/unordered_map.hpp>
#include "genfile/SNPDataSource.hpp"

namespace genfile {
	namespace impl {
		struct HLAIMPVariantDataReader ;
	}
	
	// This class represents a SNPDataSource which reads its data
	// from a HLAIMP output file.
	// Such files contain imputed genotype probabilities for a single HLA variant (which has multiple alleles).
	// This class splits these into a single biallelic variant per allele.
	// We infer allele names from the column names in the file, and set positions to dummy positions on chromosome 6.
	class HLAIMPAsBiallelicVariantDataSource: public IdentifyingDataCachingSNPDataSource
	{
	public:
		HLAIMPAsBiallelicVariantDataSource( std::string const& filename ) ;
		HLAIMPAsBiallelicVariantDataSource( std::istream& stream ) ;

	public:
		Metadata get_metadata() const ;

		unsigned int number_of_samples() const ;
		void get_sample_ids( GetSampleIds ) const ;

		OptionalSnpCount total_number_of_snps() const ;
		operator bool() const ;

		std::string get_source_spec() const ;

	private:
		std::string m_filename ;
		std::vector< std::string > m_samples ;
		std::vector< std::string > m_alleles ;
		std::vector< std::vector< std::vector< std::size_t > > > m_allele_dosage_columns ;
		std::vector< Eigen::VectorXd > m_data ;
		std::size_t m_allele_index ;
		bool m_exhausted ;

	private:
		void reset_to_start_impl() ;
		
		void read_snp_identifying_data_impl( 
			uint32_t* number_of_samples, // number_of_samples is unused.
			std::string* SNPID,
			std::string* RSID,
			Chromosome* chromosome,
			uint32_t* SNP_position,
			std::string* allele1,
			std::string* allele2
		) ;

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() ;
		void setup( std::string const& filename ) ;
		void setup( std::istream& stream ) ;
		
		friend struct impl::HLAIMPVariantDataReader ;
	} ;
}

#endif
