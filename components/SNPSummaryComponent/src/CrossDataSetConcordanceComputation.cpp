
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

	void CrossDataSetConcordanceComputation::set_comparer( genfile::SNPIdentifyingData::CompareFields const& comparer ) {
		m_comparer = comparer ;
	}

	void CrossDataSetConcordanceComputation::set_alternate_dataset(
		genfile::CohortIndividualSource::UniquePtr samples,
		std::string const& comparison_dataset_sample_id_column,
		genfile::SNPDataSource::UniquePtr snps
	) {
		m_sample_mapper.set_alternate_dataset( samples, comparison_dataset_sample_id_column ) ;
		m_alt_dataset_snps = snps ;
	}

	void CrossDataSetConcordanceComputation::operator()(
		SNPIdentifyingData const& snp,
		Genotypes const& genotypes,
		SampleSexes const&,
		genfile::VariantDataReader&,
		ResultCallback callback
	) {
		using genfile::string_utils::to_string ;
		SNPIdentifyingData alt_snp ;
		if( m_alt_dataset_snps->get_next_snp_matching( &alt_snp, snp, m_comparer )) {
			callback( "compared_variant_rsid", alt_snp.get_rsid() ) ;
			genfile::vcf::ThreshholdingGenotypeSetter< Eigen::VectorXd > setter( m_alt_dataset_genotypes, m_call_threshhold, -1, 0, 1, 2 ) ;
			m_alt_dataset_snps->read_variant_data()->get( "genotypes", setter ) ;

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
				std::string const stub = m_sample_mapper.dataset1_sample_ids()[ i->first ].as< std::string >() + "(" + to_string( i->first + 1 ) + "~" + to_string( i->second + 1 ) + ")"  ;
				if( alt_genotype != -1 ) {
					++call_count ;
					for( int g = 0; g < 3; ++g ) {
						if( genotypes( i->first, g ) > m_call_threshhold ) {
							m_main_dataset_genotype_subset( alt_dataset_sample_index ) = g ;
							m_pairwise_nonmissingness( alt_dataset_sample_index ) = 1 ;

							if( g == alt_genotype ) {
								callback( stub + ":concordance", 1 ) ;
								++concordant_call_count ;
							}
							else {
								callback( stub + ":concordance", 0 ) ;
							}
						}
					}
				} else {
					callback( stub + ":concordance", genfile::MissingValue() ) ;
				}
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

	CrossDataSetHaplotypeComparisonComputation::CrossDataSetHaplotypeComparisonComputation(
		genfile::CohortIndividualSource const& samples,
		std::string const& main_dataset_sample_id_column
	):
		m_sample_mapper( samples, main_dataset_sample_id_column ),
		m_call_threshhold( 0.9 )
	{}

	void CrossDataSetHaplotypeComparisonComputation::set_comparer( genfile::SNPIdentifyingData::CompareFields const& comparer ) {
		m_comparer = comparer ;
	}

	void CrossDataSetHaplotypeComparisonComputation::set_alternate_dataset(
		genfile::CohortIndividualSource::UniquePtr samples,
		std::string const& comparison_dataset_sample_id_column,
		genfile::SNPDataSource::UniquePtr snps
	) {
		m_sample_mapper.set_alternate_dataset( samples, comparison_dataset_sample_id_column ) ;
		m_alt_dataset_snps = snps ;
		m_relative_phase.setZero( m_sample_mapper.sample_mapping().size() ) ;
	}
	
	void CrossDataSetHaplotypeComparisonComputation::operator()(
		SNPIdentifyingData const& snp,
		Genotypes const&,
		SampleSexes const&,
		genfile::VariantDataReader& data_reader,
		ResultCallback callback
	) {
		using genfile::string_utils::to_string ;
		SNPIdentifyingData alt_snp ;
		if( m_alt_dataset_snps->get_next_snp_matching( &alt_snp, snp, m_comparer )) {
			callback( "CrossDataSetHaplotypeComparisonComputation:compared_variant_rsid", alt_snp.get_rsid() ) ;

			{
				genfile::vcf::PhasedGenotypeSetter< Eigen::MatrixXd > setter( m_haplotypes1, m_nonmissingness1 ) ;
				data_reader.get( "genotypes", setter ) ;
			}
			{
				genfile::vcf::PhasedGenotypeSetter< Eigen::MatrixXd > setter( m_haplotypes2, m_nonmissingness2 ) ;
				m_alt_dataset_snps->read_variant_data()->get( "genotypes", setter ) ;
			}
			
			// Threshhold the calls and set to missing any rows not meeting the
			// threshhold.
			m_haplotypes1.array() *= ( m_haplotypes1.array() >= m_call_threshhold ).cast< double >() ; 
			m_haplotypes2.array() *= ( m_haplotypes2.array() >= m_call_threshhold ).cast< double >() ; 

			// FIX MISSINGNESS HERE

			CrossDataSetSampleMapper::SampleMapping::const_iterator i = m_sample_mapper.sample_mapping().begin() ;
			CrossDataSetSampleMapper::SampleMapping::const_iterator const end_i = m_sample_mapper.sample_mapping().end() ;
			int call_count = 0 ;
			int concordant_call_count = 0 ;
			int het_call_count = 0 ;
			int switch_error_count = 0 ;
			for( std::size_t count = 0; i != end_i; ++i, ++count ) {
				std::string const stub = (
					m_sample_mapper.dataset1_sample_ids()[ i->first ].as< std::string >()
					+ "("
					+ to_string( i->first + 1 )
					+ "~"
					+ to_string( i->second + 1 )
					+ ")"
				) ;

				if( m_nonmissingness1.row( i->first ).sum() == 0 && m_nonmissingness2.row( i->second ).sum() == 0 ) {
					++call_count ;

					Eigen::MatrixXd::RowXpr const& h1 = m_haplotypes1.row( i->first ) ;
					Eigen::MatrixXd::RowXpr const& h2 = m_haplotypes2.row( i->second ) ;
					int const concordant = ( h1.sum() == h2.sum() ) ;
					callback( stub + ":concordance", concordant ) ;
					concordant_call_count += concordant ;
					
					int const heterozygote = ( h1.sum() == 1 ) ;
					het_call_count += heterozygote ;
					if( concordant && heterozygote ) {
						
						int const relative_phase = ( h1 == h2 ) ? 1 : -1 ;
						int const switch_error = ( m_relative_phase( count ) != 0 ) && ( m_relative_phase( count ) != relative_phase ) ;
						callback( stub + "switch_error", switch_error ) ;
						// take account of the switch for the next SNP.
						m_relative_phase( count ) = relative_phase ;
						switch_error_count += switch_error ;
					}
				} else {
					callback( stub + ":concordance", genfile::MissingValue() ) ;
					callback( stub + ":switch_error", genfile::MissingValue() ) ;
				}
			}
			callback( "pairwise non-missing haplotypes", call_count ) ;
			callback( "pairwise concordant haplotypes", concordant_call_count ) ;
			callback( "concordant heterozygous haplotypes", het_call_count ) ;
			callback( "concordant heterozygoes haplotypes with switch error", het_call_count ) ;
		}
	}

	std::string CrossDataSetHaplotypeComparisonComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "CrossDataSetHaplotypeComparisonComputation" ;
	}
}

