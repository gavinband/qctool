
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/bind.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/AlleleFlippingSNPDataSource.hpp"
#include "genfile/AlleleFlippingVariantDataReader.hpp"
#include "genfile/get_set.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	std::pair< std::vector< SNPIdentifyingData >, AlleleFlippingSNPDataSource::AlleleFlipSpec >
	AlleleFlippingSNPDataSource::get_allele_flip_spec(
		std::vector< SNPIdentifyingData > reference_snps,
		std::vector< SNPIdentifyingData > snps_to_match,
		SNPIdentifyingData::CompareFields const& comparator
	) {
		std::sort( reference_snps.begin(), reference_snps.end(), comparator ) ;
		AlleleFlipSpec allele_flips( snps_to_match.size() ) ;
		
		for( std::size_t i = 0; i < snps_to_match.size(); ++i ) {
			SNPIdentifyingData swapped_snp = snps_to_match[i] ;
			std::string const first_allele = swapped_snp.get_first_allele() ;
			swapped_snp.set_first_allele( swapped_snp.get_second_allele() ) ;
			swapped_snp.set_second_allele( first_allele ) ;

			if( snps_to_match[i].get_first_allele() == "?" || snps_to_match[i].get_second_allele() == "?" ) {
				allele_flips[i] = eUnknownFlip ;
			}
			else if(
				std::binary_search(
					reference_snps.begin(),
					reference_snps.end(),
					snps_to_match[i],
					comparator
				)
			) {
				allele_flips[i] = eNoFlip ;
			}
			else if(
				std::binary_search(
					reference_snps.begin(),
					reference_snps.end(),
					swapped_snp,
					comparator
				)
			) {
				snps_to_match[i] = swapped_snp ;
				allele_flips[i] = eFlip ;
			}
			else {
				snps_to_match[i].set_first_allele( "?" ) ;
				snps_to_match[i].set_second_allele( "?" ) ;
				allele_flips[i] = eUnknownFlip ;
			}
		}
		return std::make_pair( snps_to_match, allele_flips ) ;
	}
	
	AlleleFlippingSNPDataSource::UniquePtr AlleleFlippingSNPDataSource::create(
		SNPDataSource::UniquePtr source,
		AlleleFlipSpec const& allele_flips
	) {
		return UniquePtr( new AlleleFlippingSNPDataSource( source, allele_flips )) ;
	}
	
	AlleleFlippingSNPDataSource::AlleleFlippingSNPDataSource(
		SNPDataSource::UniquePtr source,
		AlleleFlipSpec const& allele_flips
	):
		m_source( source ),
		m_allele_flips( allele_flips )
	{
		assert( !m_source->total_number_of_snps() || m_allele_flips.size() == *m_source->total_number_of_snps() ) ;
	}
	
	SNPDataSource::Metadata AlleleFlippingSNPDataSource::get_metadata() const {
		return get_parent_source().get_metadata() ;
	}
	
	std::string AlleleFlippingSNPDataSource::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return m_source->get_summary( prefix + "  " ) ;
	}

	void AlleleFlippingSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		char allele_flip = eUnknownFlip ;
		if( number_of_snps_read() < m_allele_flips.size() ) {
			allele_flip = m_allele_flips[ number_of_snps_read() ] ;
		}
		switch( allele_flip ) {
			case eNoFlip:
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
			case eFlip:
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
			case eUnknownFlip:
				m_source->get_snp_identifying_data(
					set_number_of_samples,
					set_SNPID,
					set_RSID,
					set_chromosome,
					set_SNP_position,
					ignore(),
					ignore()
				) ;
				set_allele1( "?" ) ;
				set_allele2( "?" ) ;
				break ;
			default:
				assert(0) ;
		}
	}

	VariantDataReader::UniquePtr AlleleFlippingSNPDataSource::read_variant_data_impl() {
		return VariantDataReader::UniquePtr(
			new AlleleFlippingVariantDataReader(
				number_of_samples(),
				m_source->read_variant_data(),
				m_allele_flips[ m_source->number_of_snps_read() ]
			)
		) ;
	}
	
	AlleleFlippingSNPDataSource::AlleleFlippingGenotypeProbabilitySetter::AlleleFlippingGenotypeProbabilitySetter( GenotypeProbabilitySetter setter ):
		m_setter( setter )
	{}

	void AlleleFlippingSNPDataSource::AlleleFlippingGenotypeProbabilitySetter::operator()( std::size_t i, double AA, double AB, double BB ) const {
		m_setter( i, BB, AB, AA ) ;
	}

	void AlleleFlippingSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}

	void AlleleFlippingSNPDataSource::reset_to_start_impl() {
		m_source->reset_to_start() ;
	}
}
