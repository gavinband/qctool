#ifndef GENFILE_STRAND_ALIGNING_SNP_DATA_SOURCE_HPP
#define GENFILE_STRAND_ALIGNING_SNP_DATA_SOURCE_HPP

#include <vector>
#include <map>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"

namespace genfile {
	struct StrandAligningSNPDataSource: public SNPDataSource
	{
	public:
		typedef std::vector< char > StrandAlignments ;
		static char const eUnknownStrand = '?' ;
		static char const eForwardStrand = '+' ;
		static char const eReverseStrand  = '-' ;

		static std::pair< std::vector< SNPIdentifyingData >, StrandAlignments > create_strand_alignments(
			std::vector< SNPIdentifyingData > snps,
			std::map< SNPIdentifyingData, char > known_strand_alignments
		) ;
		
	public:
		typedef std::auto_ptr< StrandAligningSNPDataSource > UniquePtr ;
		static UniquePtr create( SNPDataSource::UniquePtr source, StrandAlignments const& strand_alignments ) ;

		StrandAligningSNPDataSource( SNPDataSource::UniquePtr source, StrandAlignments const& strand_alignments ) ;

		std::vector< SNPIdentifyingData > const& get_aligned_snps() const { return m_aligned_snps ; }
	private:
		
		SNPDataSource::UniquePtr m_source ;
		StrandAlignments const m_strand_alignments ;
		std::vector< SNPIdentifyingData > const m_aligned_snps ;
		
	public:
		operator bool() const { return *m_source ; }
		unsigned int number_of_samples() const { return m_source->number_of_samples() ; }
		unsigned int total_number_of_snps() const { return m_source->total_number_of_snps() ; }
		std::string get_source_spec() const { return "strand-aligned:" + m_source->get_source_spec() ; }
		SNPDataSource const& get_parent_source() const {
			return *m_source ;
		}
		SNPDataSource const& get_base_source() const {
			return m_source->get_base_source() ;
		}

		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;

	private:

		struct StrandFlippingGenotypeProbabilitySetter
		{
			StrandFlippingGenotypeProbabilitySetter( GenotypeProbabilitySetter setter ) ;
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

		void read_snp_probability_data_impl(
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) ;

		void ignore_snp_probability_data_impl() ;
		void reset_to_start_impl() ;
	} ;
}

#endif
