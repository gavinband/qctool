#include <vector>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/StrandAligningSNPDataSource.hpp"
#include "genfile/get_set.hpp"

namespace genfile {
	std::pair< std::vector< SNPIdentifyingData >, StrandAligningSNPDataSource::StrandAlignments > StrandAligningSNPDataSource::create_strand_alignments(
		std::vector< SNPIdentifyingData > snps,
		std::map< SNPIdentifyingData, char > known_strand_alignments
	) {
		StrandAlignments result( snps.size() ) ;
		for( std::size_t snp_i = 0; snp_i < snps.size(); ++snp_i) {
			std::map< SNPIdentifyingData, char >::const_iterator where = known_strand_alignments.find( snps[ snp_i ] ) ;
			if( where == known_strand_alignments.end() ) {
				result[ snp_i ] = eUnknownStrand ;
				snps[ snp_i ].first_allele() = snps[ snp_i ].second_allele() = '?' ;
			}
			else {
				switch( where->second ) {
					case eForwardStrand:
						result[ snp_i ] = eForwardStrand ;
						break ;
					case eReverseStrand:
						result[ snp_i ] = eReverseStrand ;
						std::swap( snps[ snp_i ].first_allele(), snps[ snp_i ].second_allele() ) ;
						break ;
					case eUnknownStrand:
						result[ snp_i ] = eUnknownStrand ;
						snps[ snp_i ].first_allele() = '?' ;
						snps[ snp_i ].second_allele() = '?' ;
						break ;
					default:
						assert(0);
						break ;
				}
			}
		}
		return std::make_pair( snps, result ) ;
	}
	
	StrandAligningSNPDataSource::UniquePtr StrandAligningSNPDataSource::create(
		SNPDataSource::UniquePtr source,
		StrandAlignments const& strand_alignments
	) {
		return UniquePtr( new StrandAligningSNPDataSource( source, strand_alignments )) ;
	}
	
	StrandAligningSNPDataSource::StrandAligningSNPDataSource(
		SNPDataSource::UniquePtr source,
		StrandAlignments const& strand_alignments
	):
		m_source( source ),
		m_strand_alignments( strand_alignments )
	{
		assert( m_strand_alignments.size() == m_source->total_number_of_snps() ) ;
	}
	
	std::string StrandAligningSNPDataSource::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return m_source->get_summary( prefix + "  " ) ;
	}

	void StrandAligningSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		char strand_alignment = eUnknownStrand ;
		if( number_of_snps_read() < m_strand_alignments.size() ) {
			strand_alignment = m_strand_alignments[ number_of_snps_read() ] ;
		}
		switch( strand_alignment ) {
			case eForwardStrand:
				m_source->get_snp_identifying_data(
					set_number_of_samples,
					set_SNPID,
					set_RSID,
					set_chromosome,
					set_SNP_position,
					set_allele1,
					set_allele2
				) ;
				break ;
			case eReverseStrand:
				m_source->get_snp_identifying_data(
					set_number_of_samples,
					set_SNPID,
					set_RSID,
					set_chromosome,
					set_SNP_position,
					set_allele2,
					set_allele1
				) ;
				break ;
			case eUnknownStrand:
				m_source->get_snp_identifying_data(
					set_number_of_samples,
					set_SNPID,
					set_RSID,
					set_chromosome,
					set_SNP_position,
					ignore(),
					ignore()
				) ;
				set_allele1( '?' ) ;
				set_allele2( '?' ) ;
				break ;
			default:
				assert(0) ;
		}
	}

	void StrandAligningSNPDataSource::read_snp_probability_data_impl(
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		switch( m_strand_alignments[ m_source->number_of_snps_read() ] ) {
			case eForwardStrand:
				m_source->read_snp_probability_data( set_genotype_probabilities ) ;
				break ;
			case eReverseStrand:
				m_source->read_snp_probability_data(
					StrandFlippingGenotypeProbabilitySetter( set_genotype_probabilities )
				);
				break ;
			case eUnknownStrand:
				m_source->ignore_snp_probability_data() ;
				for( std::size_t i = 0; i < number_of_samples(); ++i ) {
					set_genotype_probabilities( i, 0.0, 0.0, 0.0 ) ;
				}
				break ; 
			default:
				assert(0) ; 
		}
	}

	StrandAligningSNPDataSource::StrandFlippingGenotypeProbabilitySetter::StrandFlippingGenotypeProbabilitySetter( GenotypeProbabilitySetter setter ):
		m_setter( setter )
	{}

	void StrandAligningSNPDataSource::StrandFlippingGenotypeProbabilitySetter::operator()( std::size_t i, double AA, double AB, double BB ) const {
		m_setter( i, BB, AB, AA ) ;
	}

	void StrandAligningSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}

	void StrandAligningSNPDataSource::reset_to_start_impl() {
		m_source->reset_to_start() ;
	}
}
