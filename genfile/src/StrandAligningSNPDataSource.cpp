
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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
				default:
					// Anything else we'll just return verbatim, this handles indels and deletions (I and D).
					// But it would not handle complex alleles made up of base sequences.
					return allele ;
					break ;
			}
		}
		
		std::string complement( std::string const& allele ) {
			if( allele.size() != 1 ) {
				throw BadArgumentError( "StrandAligningSNPDataSource::complement", "allele=" + allele ) ;
			}
			return std::string( 1, complement( allele[0] ) ) ;
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
				snps[ snp_i ].first_allele() = ( snps[ snp_i ].get_first_allele() == impl::complement( snps[ snp_i ].get_first_allele() ) ) ? snps[ snp_i ].get_first_allele() : "?" ;
				snps[ snp_i ].second_allele() = ( snps[ snp_i ].get_second_allele() == impl::complement( snps[ snp_i ].get_second_allele() ) ) ? snps[ snp_i ].get_second_allele() : "?" ;
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
						snps[ snp_i ].first_allele() = ( snps[ snp_i ].get_first_allele() == impl::complement( snps[ snp_i ].get_first_allele() ) ) ? snps[ snp_i ].get_first_allele() : "?" ;
						snps[ snp_i ].second_allele() = ( snps[ snp_i ].get_second_allele() == impl::complement( snps[ snp_i ].get_second_allele() ) ) ? snps[ snp_i ].get_second_allele() : "?" ;
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
		std::string allele1, allele2 ;
		
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
					set_value( allele1 ),
					set_value( allele2 )
				) ;
				if( allele1 == impl::complement( allele1 ) ) {
					set_allele1( allele1 ) ;
				}
				else {
					set_allele1( "?" ) ;
				}

				if( allele2 == impl::complement( allele2 ) ) {
					set_allele2( allele2 ) ;
				}
				else {
					set_allele2( "?" ) ;
				}
				break ;
			default:
				assert(0) ;
		}
	}

 	VariantDataReader::UniquePtr StrandAligningSNPDataSource::read_variant_data_impl() {
		return m_source->read_variant_data() ;
	}

	void StrandAligningSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}

	void StrandAligningSNPDataSource::reset_to_start_impl() {
		m_source->reset_to_start() ;
	}
}
