
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_MERGING_SNP_DATA_SOURCE_HPP
#define GENFILE_MERGING_SNP_DATA_SOURCE_HPP

#include <vector>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPDataSource.hpp"

namespace genfile {
	struct MergingSNPDataSource: public SNPDataSource {

		typedef std::auto_ptr< MergingSNPDataSource > UniquePtr ;
		static std::vector< std::string > get_merge_strategies() ;
		static UniquePtr create(
			std::string const& merge_strategy,
			SNPIdentifyingData::CompareFields const& = SNPIdentifyingData::CompareFields()
		) ;

		MergingSNPDataSource( genfile::SNPIdentifyingData::CompareFields const& ) ;
		virtual ~MergingSNPDataSource() ;
		
		void add_source( SNPDataSource::UniquePtr, std::string const& id_prefix = "" ) ;

		operator bool() const ;
		// Return the number of samples represented in the snps in this source.
		unsigned int number_of_samples() const ;
		// Return the total number of snps the source contains.
		OptionalSnpCount total_number_of_snps() const ;

		// Return a string identifying the source of the SNP data
		std::string get_source_spec() const ;

	private:
		void get_snp_identifying_data_impl( 
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			ChromosomeSetter const& set_chromosome,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) ;	

		VariantDataReader::UniquePtr read_variant_data_impl()  ;

		void ignore_snp_probability_data_impl() ;
		void reset_to_start_impl() ;
	
	private:
		std::vector< SNPDataSource* > m_sources ;
		typedef std::multimap< SNPIdentifyingData, std::size_t, SNPIdentifyingData::CompareFields > CurrentSnps ;
		CurrentSnps m_current_snps ;
		std::vector< std::string > m_merge_id_prefixes ;
		void get_top_snp_in_source( std::size_t source_i ) ;
	protected:
		CurrentSnps& current_snps() { return m_current_snps ; }
		SNPDataSource& get_source( std::size_t i ) { return *m_sources[i] ; }
		void discard_top_snp_and_get_next_candidate() ;
		virtual void discard_top_snp_and_get_next() = 0 ;
	} ;
	
	struct DropDuplicatesStrategyMergingSNPDataSource: public MergingSNPDataSource {
		DropDuplicatesStrategyMergingSNPDataSource( SNPIdentifyingData::CompareFields const& compare_fields ) ;
		void discard_top_snp_and_get_next() ;
	} ;

	struct KeepAllStrategyMergingSNPDataSource: public MergingSNPDataSource {
		KeepAllStrategyMergingSNPDataSource( SNPIdentifyingData::CompareFields const& compare_fields ) ;
		void discard_top_snp_and_get_next() ;
	} ;
}

#endif
