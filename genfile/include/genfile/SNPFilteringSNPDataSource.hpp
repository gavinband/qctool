
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SNPFILTERINGSNPDATASOURCE_HPP
#define GENFILE_SNPFILTERINGSNPDATASOURCE_HPP

#include <vector>
#include <set>
#include <memory>
#include "genfile/Error.hpp"
#include "genfile/SNPDataSource.hpp"

namespace genfile {
	// class SNPFilteringSNPDataSource
	// This class represents a view of its parent SNPDataSource such that
	// only SNPs whose chromosome, position, id fields and alleles satisfy the
	// given test are visible.
	class SNPFilteringSNPDataSource: public SNPDataSource
	{
	public:
		typedef std::auto_ptr< SNPFilteringSNPDataSource > UniquePtr ;
		typedef std::vector< std::size_t > IndexList ;
	public:
		// Create a SNPFilteringSNPDataSource from the given source and the given sample indices.
		static UniquePtr create( SNPDataSource::UniquePtr source, IndexList const& indices_of_snps_to_include ) ;
		
	public:
		SNPFilteringSNPDataSource( SNPDataSource::UniquePtr source, IndexList indices_of_snps_to_include ) ;

		operator bool() const ;
		Metadata get_metadata() const ;
		unsigned int number_of_samples() const ;
		void get_sample_ids( GetSampleIds ) const ;
		OptionalSnpCount total_number_of_snps() const ;
		OptionalSnpCount total_number_of_snps_before_filtering() const ;
		std::string get_source_spec() const ;

		SNPDataSource const& get_parent_source() const ;
		SNPDataSource const& get_base_source() const ;
		
	private:
		void reset_to_start_impl() ;

		void get_snp_identifying_data_impl( 
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			ChromosomeSetter const& set_chromosome,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) ;

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() ;

	protected:
		SNPDataSource& source() { return *m_source ; }
		
	private:
		
		SNPDataSource::UniquePtr m_source ;
		std::set< std::size_t > m_indices_of_excluded_snps ;
	} ;
}

#endif
