
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>
#include <algorithm>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/StrandAligningSNPDataSource.hpp"
#include "genfile/OffsetFlippedAlleleSetter.hpp"
#include "genfile/AlleleFlippingVariantDataReader.hpp"
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
			std::string result = allele ;
			for( std::size_t i = 0; i < result.size(); ++i ) {
				result[i] = complement( result[i] ) ;
			}
			return result ;
		}
	}
	
	std::string StrandAligningSNPDataSource::apply_strand( std::string const& allele, char strand ) {
		assert( strand == StrandAligningSNPDataSource::eForwardStrand || strand == StrandAligningSNPDataSource::eReverseStrand ) ;
		if( strand == StrandAligningSNPDataSource::eForwardStrand ) {
			return allele ;
		} else {
			return impl::complement( allele ) ;
		}
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
		m_strand_alignments( strand_alignments ),
		m_include_unknown_strand_or_flip( false )
	{}
	
	SNPDataSource::Metadata StrandAligningSNPDataSource::get_metadata() const {
		return m_source->get_metadata() ;
	}
	
	std::string StrandAligningSNPDataSource::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return m_source->get_summary( prefix + "  " ) ;
	}

	void StrandAligningSNPDataSource::get_snp_identifying_data_impl( 
VariantIdentifyingData* variant	) {
		SNPIdentifyingData source_snp ;
		for( ; m_source->get_snp_identifying_data( source_snp ); m_source->ignore_snp_probability_data() ) {
			m_current_strand_flip_spec = get_strand_alignment( source_snp ) ;
			//std::cerr << "looking at SNP " << source_snp << ".\n" ;
			if( 
				m_include_unknown_strand_or_flip ||
				( m_current_strand_flip_spec.strand != eUnknownStrand && m_current_strand_flip_spec.flip != eUnknownFlip )
			) {
				char strand_alignment = m_current_strand_flip_spec.strand ;

				std::string allele1, allele2 ;

				m_source->get_snp_identifying_data(
					set_number_of_samples,
					set_SNPID,
					set_RSID,
					set_chromosome,
					set_SNP_position,
					set_value( allele1 ),
					set_value( allele2 )
				) ;

				switch( strand_alignment ) {
					case eForwardStrand:
						//ok.
						break ;
					case eReverseStrand:
						//ok.
						allele1 = impl::complement( allele1 ) ;
						allele2 = impl::complement( allele2 ) ;
						break ;
					case eUnknownStrand:
						if( allele1 != impl::complement( allele1 ) ) {
							allele1 += ( "/" + impl::complement( allele1 ) ) ;
						}
						if( allele2 != impl::complement( allele2 ) ) {
							allele2 += ( "/" + impl::complement( allele2 ) ) ;
						}
						break ;
					default:
						assert(0) ;
				}

				switch( m_current_strand_flip_spec.flip ) {
					case eNoFlip:
						// no change needed
						break ;
					case eFlip:
						std::swap( allele1, allele2 ) ;
						break ;
					case eUnknownFlip:
						allele1 = allele1 + "/" + allele2 ;
						allele2 = allele2 + "/" + source_snp.get_first_allele() ;
						break ;
					default:
						assert(0) ;
				}
				
				set_allele1( allele1 ) ;
				set_allele2( allele2 ) ;
				return ;
			}
		}
	}

	StrandAligningSNPDataSource::StrandFlipSpec StrandAligningSNPDataSource::get_strand_alignment( SNPIdentifyingData const& snp ) const {
		StrandFlipSpec result = StrandFlipSpec() ;
		StrandAlignments::const_iterator where = m_strand_alignments.find( snp ) ;
		if( where != m_strand_alignments.end() ) {
			result = where->second ;
		}
		return result ;
	}

 	VariantDataReader::UniquePtr StrandAligningSNPDataSource::read_variant_data_impl() {
		return VariantDataReader::UniquePtr(
			new AlleleFlippingVariantDataReader(
				number_of_samples(),
				m_source->read_variant_data(),
				m_current_strand_flip_spec.flip
			)
		) ;
	}

	void StrandAligningSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}

	void StrandAligningSNPDataSource::reset_to_start_impl() {
		m_source->reset_to_start() ;
	}
}
