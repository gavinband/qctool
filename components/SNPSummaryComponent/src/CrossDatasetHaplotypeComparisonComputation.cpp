
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
#include "components/SNPSummaryComponent/CrossDataSetHaplotypeComparisonComputation.hpp"

namespace snp_stats {
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

	void CrossDataSetHaplotypeComparisonComputation::set_match_alleles() {
		m_match_alleles = true ;
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
			{
				genfile::vcf::PhasedGenotypeSetter< Eigen::MatrixXd > setter( m_haplotypes1, m_nonmissingness1 ) ;
				data_reader.get( "genotypes", setter ) ;
			}
			{
				double A_coding = 0 ;
				double B_coding = 1 ;
				bool matching_alleles = true ;
				if( m_match_alleles ) {
					if( alt_snp.get_first_allele() == snp.get_first_allele() && alt_snp.get_second_allele() == snp.get_second_allele() ) {
						A_coding = 0 ;
						B_coding = 1 ;
					}
					else if( alt_snp.get_first_allele() == snp.get_second_allele() && alt_snp.get_second_allele() == snp.get_first_allele() ) {
						A_coding = 1 ;
						B_coding = 0 ;
					} else {
						matching_alleles = false ;
					}
				}
				genfile::vcf::PhasedGenotypeSetter< Eigen::MatrixXd > setter( m_haplotypes2, m_nonmissingness2, 0, A_coding, B_coding ) ;
				m_alt_dataset_snps->read_variant_data()->get( "genotypes", setter ) ;

				if( !matching_alleles ) {
					callback( "comment", "Alleles in main and comparison datasets do not match." ) ;
					// alleles don't match, don't  bother doing any computation.
					return ;
				}
			}
			
			callback( "CrossDataSetHaplotypeComparisonComputation:compared_variant_rsid", alt_snp.get_rsid() ) ;
			callback( "compared_variant_alleleA", alt_snp.get_first_allele() ) ;
			callback( "compared_variant_alleleB", alt_snp.get_second_allele() ) ;

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

				if( m_nonmissingness1.row( i->first ).sum() == 2 && m_nonmissingness2.row( i->second ).sum() == 2 ) {
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
					} else {
						callback( stub + "switch_error", genfile::MissingValue() ) ;
					}
				} else {
					callback( stub + ":concordance", genfile::MissingValue() ) ;
					callback( stub + ":switch_error", genfile::MissingValue() ) ;
				}
			}
			callback( "pairwise_non-missing_haplotypes", call_count ) ;
			callback( "pairwise_concordant_haplotypes", concordant_call_count ) ;
			callback( "concordant_heterozygous_haplotypes", het_call_count ) ;
			callback( "concordant_heterozygous_haplotypes_with_switch_error", switch_error_count ) ;
		}
	}

	std::string CrossDataSetHaplotypeComparisonComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "CrossDataSetHaplotypeComparisonComputation" ;
	}
}

