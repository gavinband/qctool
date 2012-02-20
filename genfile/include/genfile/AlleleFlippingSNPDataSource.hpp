#ifndef GENFILE_ALLELE_FLIPPING_SNP_DATA_SOURCE_HPP
#define GENFILE_ALLELE_FLIPPING_SNP_DATA_SOURCE_HPP

#include <vector>
#include <map>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"

namespace genfile {
	namespace impl {
		struct AlleleFlippingSNPDataReader ;
	}
	struct AlleleFlippingSNPDataSource: public SNPDataSource
	{
		friend class impl::AlleleFlippingSNPDataReader ;
	public:
		typedef std::vector< char > AlleleFlipSpec ;
		static char const eUnknownFlip = '?' ;
		static char const eNoFlip = '+' ;
		static char const eFlip  = '-' ;

		static std::pair< std::vector< SNPIdentifyingData >, AlleleFlipSpec > get_allele_flip_spec(
			std::vector< SNPIdentifyingData > reference_snps,
			std::vector< SNPIdentifyingData > snps_to_match,
			SNPIdentifyingData::CompareFields const&
		) ;
		
	public:
		typedef std::auto_ptr< AlleleFlippingSNPDataSource > UniquePtr ;
		static UniquePtr create( SNPDataSource::UniquePtr source, AlleleFlipSpec const& allele_flips ) ;
		
		AlleleFlippingSNPDataSource( SNPDataSource::UniquePtr source, AlleleFlipSpec const& allele_flips ) ;
		
		std::vector< SNPIdentifyingData > const& get_flipped_snps() const { return m_flipped_snps ; }
	private:
		
		SNPDataSource::UniquePtr m_source ;
		AlleleFlipSpec const m_allele_flips ;
		std::vector< SNPIdentifyingData > const m_flipped_snps ;
		
	public:
		operator bool() const { return *m_source ; }
		unsigned int number_of_samples() const { return m_source->number_of_samples() ; }
		OptionalSnpCount total_number_of_snps() const { return m_source->total_number_of_snps() ; }
		std::string get_source_spec() const { return "allele-flipped:" + m_source->get_source_spec() ; }
		SNPDataSource const& get_parent_source() const {
			return *m_source ;
		}
		SNPDataSource const& get_base_source() const {
			return m_source->get_base_source() ;
		}

		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;

	private:

		struct AlleleFlippingGenotypeProbabilitySetter
		{
			AlleleFlippingGenotypeProbabilitySetter( GenotypeProbabilitySetter setter ) ;
			void operator()( std::size_t i, double AA, double AB, double BB ) const ;
		private:
			GenotypeProbabilitySetter m_setter ;
		} ;

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
		void reset_to_start_impl() ;
	} ;
}

#endif
