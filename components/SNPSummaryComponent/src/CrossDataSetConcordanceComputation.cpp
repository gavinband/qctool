
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <boost/bind.hpp>
#include <Eigen/Core>
#include "genfile/VariantEntry.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSource.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "components/SNPSummaryComponent/CrossDataSetConcordanceComputation.hpp"
#include "metro/correlation.hpp"

namespace snp_stats {
	namespace impl {
		CrossDataSetSampleMapper::SampleMapping create_sample_mapping(
			CrossDataSetSampleMapper::SampleIdList const& set1,
			CrossDataSetSampleMapper::SampleIdList const& set2
		) {
			CrossDataSetSampleMapper::SampleMapping result ;
			for( std::size_t i = 0; i < set1.size(); ++i ) {
				CrossDataSetSampleMapper::SampleIdList::const_iterator where = std::find( set2.begin(), set2.end(), set1[i] ) ;
				for( ; where != set2.end(); where = std::find( where, set2.end(), set1[i] ) ) {
					result.insert( std::make_pair( i, std::size_t( where - set2.begin() ) ) ) ;
					++where ;
				}
			}
			return result ;
		}
	}

	CrossDataSetSampleMapper::CrossDataSetSampleMapper(
		genfile::CohortIndividualSource const& dataset1_samples,
		std::string const& sample_id_column
	):
		m_dataset1_samples( dataset1_samples )
	{
		m_dataset1_samples.get_column_values( "ID_1", boost::bind( &SampleIdList::push_back, &m_dataset1_main_ids, _2 ) ) ;
		m_dataset1_samples.get_column_values( sample_id_column, boost::bind( &SampleIdList::push_back, &m_dataset1_sample_ids, _2 ) ) ;
	}

	void CrossDataSetSampleMapper::set_alternate_dataset(
		genfile::CohortIndividualSource::UniquePtr samples,
		std::string const& sample_id_column
	) {
		m_dataset2_samples = samples ;
		m_dataset2_samples->get_column_values( sample_id_column, boost::bind( &SampleIdList::push_back, &m_dataset2_sample_ids, _2 ) ) ;
		m_sample_mapping = impl::create_sample_mapping( m_dataset1_sample_ids, m_dataset2_sample_ids ) ;
	}

	CrossDataSetConcordanceComputation::CrossDataSetConcordanceComputation(
		genfile::CohortIndividualSource const& samples,
		std::string const& main_dataset_sample_id_column
	):
		m_sample_mapper( samples, main_dataset_sample_id_column ),
		m_call_threshhold( 0.9 )
	{}

	void CrossDataSetConcordanceComputation::set_comparer( genfile::VariantIdentifyingData::CompareFields const& comparer ) {
		m_comparer = comparer ;
	}

	void CrossDataSetConcordanceComputation::set_match_alleles() {
		m_match_alleles = true ;
	}

	void CrossDataSetConcordanceComputation::set_alternate_dataset(
		genfile::CohortIndividualSource::UniquePtr samples,
		std::string const& comparison_dataset_sample_id_column,
		genfile::SNPDataSource::UniquePtr snps
	) {
		m_sample_mapper.set_alternate_dataset( samples, comparison_dataset_sample_id_column ) ;
		m_alt_dataset_snps = snps ;
	}

	namespace {
		std::string format_genotype( int genotype ) {
			if( genotype == -1 ) {
				return "NA" ;
			} else {
				return genfile::string_utils::to_string( genotype ) ;
			}
		}

	}

	void CrossDataSetConcordanceComputation::operator()(
		VariantIdentifyingData const& snp,
		Genotypes const& genotypes,
		SampleSexes const&,
		genfile::VariantDataReader&,
		ResultCallback callback
	) {
		assert( snp.number_of_alleles() == 2 ) ;
		
		using genfile::string_utils::to_string ;
		genfile::VariantIdentifyingData alt_snp ;
		if( m_alt_dataset_snps->get_next_snp_matching( &alt_snp, snp, m_comparer )) {
			// Get comparison dataset genotype data.
			{
				double AA = 0, AB = 1, BB = 2 ;
				// If necessary, match alleles to the main dataset.
				if( m_match_alleles ) {
					if( alt_snp.get_allele(0) == snp.get_allele(0) && alt_snp.get_allele(1) == snp.get_allele(1) ) {
						AA = 0 ;
						AB = 1 ;
						BB = 2 ;
						// ok, coding is fine.
					}
					else if( alt_snp.get_allele(0) == snp.get_allele(1) && alt_snp.get_allele(1) == snp.get_allele(0) ) {
						AA = 2 ;
						AB = 1 ;
						BB = 0 ;
						alt_snp.swap_alleles() ;
					} else {
						// Don't know which way round alleles are.
						// Set all data to missing.
						AA = -1 ;
						AB = -1 ;
						BB = -1 ;
					}
				}
				m_alt_dataset_snps->read_variant_data()->get(
					":genotypes:",
					genfile::vcf::get_threshholded_calls( m_alt_dataset_genotypes, m_call_threshhold, -1, AA, AB, BB )
				) ;
			}
			
			callback( "compared_variant_rsid", alt_snp.get_primary_id() ) ;
			callback( "compared_variant_alleleA", alt_snp.get_allele(0) ) ;
			callback( "compared_variant_alleleB", alt_snp.get_allele(1) ) ;
			
			CrossDataSetSampleMapper::SampleMapping::const_iterator i = m_sample_mapper.sample_mapping().begin() ;
			CrossDataSetSampleMapper::SampleMapping::const_iterator const end_i = m_sample_mapper.sample_mapping().end() ;

			int call_count = 0 ;
			int concordant_call_count = 0 ;
			
			m_main_dataset_genotype_subset.resize( m_alt_dataset_genotypes.size() ) ;
			m_pairwise_nonmissingness.resize( m_alt_dataset_genotypes.size() ) ;
			m_pairwise_nonmissingness.setZero() ;

			for( std::size_t count = 0; i != end_i; ++i, ++count ) {
				std::size_t const alt_dataset_sample_index = i->second ;
				double const alt_genotype = m_alt_dataset_genotypes( alt_dataset_sample_index ) ;
				double main_genotype = std::numeric_limits< double >::quiet_NaN() ;
				if( genotypes.row( i->first ).maxCoeff( &main_genotype ) < m_call_threshhold ) {
					main_genotype = -1 ;
					m_pairwise_nonmissingness( alt_dataset_sample_index ) = 0 ;
				} else {
					m_main_dataset_genotype_subset( alt_dataset_sample_index ) = main_genotype ;
					m_pairwise_nonmissingness( alt_dataset_sample_index ) = 1 ;
				}
				std::string const stub = m_sample_mapper.dataset1_main_ids()[ i->first ].as< std::string >() + "(" + to_string( i->first + 1 ) + "~" + to_string( i->second + 1 ) + ")"  ;
				if( main_genotype != -1 && alt_genotype != -1 ) {
					++call_count ;
					concordant_call_count += ( main_genotype == alt_genotype ) ? 1 : 0 ;
				}
				callback( stub + ":genotype", format_genotype( main_genotype ) + "," + format_genotype( alt_genotype ) ) ;
			}

			callback( "pairwise non-missing calls", call_count ) ;
			callback( "pairwise concordant calls", concordant_call_count ) ;
			callback( "concordance", double( concordant_call_count ) / double( call_count ) ) ;
			callback( "correlation", metro::compute_correlation( m_main_dataset_genotype_subset, m_alt_dataset_genotypes, m_pairwise_nonmissingness )) ;
		}
	}

	std::string CrossDataSetConcordanceComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "CrossDataSetConcordanceComputation" ;
	}
}

