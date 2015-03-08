
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPDATASOURCERACK_HPP
#define SNPDATASOURCERACK_HPP

#include <vector>
#include <memory>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/GenomePosition.hpp"

namespace genfile {
	
	namespace impl {
		struct RackVariantDataReader ;
	}

	struct SNPDataSourceRack: public SNPDataSource
	{
	public:
		struct Error: public SNPDataError
		{
			Error( std::size_t source_index, SNPIdentifyingData snp )
				: m_source_index( source_index ),
				  m_snp( snp )
			{}

			~Error() throw() {}

			std::size_t source_index() const { return m_source_index ; }
			SNPIdentifyingData snp() const { return m_snp ; }

		private:
			std::size_t m_source_index ;
			SNPIdentifyingData m_snp ;
		} ;

		struct MissingSNPError: public Error
		{
			MissingSNPError( std::size_t source_index, SNPIdentifyingData snp )
				: Error( source_index, snp )
			{}
			~MissingSNPError() throw() {}
			char const* what() const throw() { return "MissingSNPError" ; }
			SNPIdentifyingData const& snp() const { return m_snp ; }
		private:
			SNPIdentifyingData const m_snp ;
		} ;

		struct SNPMismatchError: public Error
		{
			SNPMismatchError( std::size_t source_index, SNPIdentifyingData snp )
				: Error( source_index, snp )
			{}
			~SNPMismatchError() throw() {}
			char const* what() const throw() { return "SNPMismatchError" ; }
		} ;
		
	public:
		typedef std::auto_ptr< SNPDataSourceRack > UniquePtr ;
		static UniquePtr create(
			std::vector< wildcard::FilenameMatch > const& filenames
		) ;

		static UniquePtr create(
			std::vector< wildcard::FilenameMatch > const& filenames,
			std::string const& snp_match_fields
		) ;

	public:
		SNPDataSourceRack() ;
		SNPDataSourceRack( std::string const& snp_match_fields ) ;
		SNPDataSourceRack( SNPIdentifyingData::CompareFields const& comparator ) ;
		~SNPDataSourceRack() ;

		void add_source( std::auto_ptr< SNPDataSource > source ) ;
		SNPDataSource& get_source( std::size_t ) const ;

		Metadata get_metadata() const ;
		unsigned int number_of_samples() const ;
		OptionalSnpCount total_number_of_snps() const ;
		operator bool() const ;
		std::string get_source_spec() const ;
		std::string get_summary( std::string const& prefix, std::size_t width ) const ;

	protected:

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
	
	private:
		
		bool move_source_to_snp_matching(
			std::size_t source_i,
			SNPIdentifyingData const& reference_snp,
			SNPIdentifyingData* result
		) ;
		
		std::vector< SNPIdentifyingData > get_intersected_snps(
			std::vector< SNPIdentifyingData > const& snps1,
			std::vector< SNPIdentifyingData > const& snps2
		) const ;

		void check_snps_are_sorted_by_position( std::vector< SNPIdentifyingData > const& snps, std::size_t cohort_index ) ;

		char const get_flip( std::size_t i ) const { return m_flips[i]; }

		struct RackGenotypeProbabilitySetter
		{
			RackGenotypeProbabilitySetter( GenotypeProbabilitySetter const& base_setter, uint32_t index_of_first_sample ) ;

			void operator()( std::size_t, double, double, double ) const ;
			
		private:
			GenotypeProbabilitySetter const& m_base_setter ;
			uint32_t const m_index_of_first_sample ;
		} ;
		
		friend struct impl::RackVariantDataReader ;

		std::vector< SNPDataSource* > m_sources ;
		std::vector< char > m_flips ;
		uint32_t m_number_of_samples ;
		bool m_read_past_end ;
		SNPIdentifyingData::CompareFields m_comparator ;
		Metadata m_metadata ;
	} ;
}

#endif
