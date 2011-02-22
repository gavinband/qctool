#include <vector>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/StrandAligningSNPDataSource.hpp"
#include "genfile/get_set.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	namespace impl {
		char complement( char allele ) {
			switch( allele ) {
				case 'A': return 'T' ; break ;
				case 'T': return 'A' ; break ;
				case 'C': return 'G' ; break ;
				case 'G': return 'C' ; break ;
				case '?': return '?' ; break ;
				case 'I':
				case 'D':
					// Some platforms have indels.  I guess these stay the same when complemented: if DNA
					// is deleted on one strand, it must be on the other as well.
					return allele ;
					break ;
				default:
					throw BadArgumentError( "StrandAligningSNPDataSource::complement", "allele=" + std::string( 1, allele ) ) ;
					break ;
			}
		}
	}

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
						snps[ snp_i ].first_allele() = impl::complement( snps[ snp_i ].get_first_allele() ) ;
						snps[ snp_i ].second_allele() = impl::complement( snps[ snp_i ].get_second_allele() ) ;
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
		char allele1, allele2 ;
		
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
					set_value( allele1 ),
					set_value( allele2 )
				) ;
				set_allele1( impl::complement( allele1 )) ;
				set_allele2( impl::complement( allele2 )) ;
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
		m_source->read_snp_probability_data( set_genotype_probabilities ) ;
	}

	void StrandAligningSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}

	void StrandAligningSNPDataSource::reset_to_start_impl() {
		m_source->reset_to_start() ;
	}
}
