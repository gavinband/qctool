
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

namespace snp_stats {
	namespace impl {
		CrossDataSetConcordanceComputation::SampleMapping create_sample_mapping(
			CrossDataSetConcordanceComputation::SampleIdList const& set1,
			CrossDataSetConcordanceComputation::SampleIdList const& set2
		) {
			CrossDataSetConcordanceComputation::SampleMapping result ;
			for( std::size_t i = 0; i < set1.size(); ++i ) {
				CrossDataSetConcordanceComputation::SampleIdList::const_iterator where = std::find( set2.begin(), set2.end(), set1[i] ) ;
				for( ; where != set2.end(); where = std::find( where, set2.end(), set1[i] ) ) {
					result.insert( std::make_pair( i, std::size_t( where - set2.begin() ) ) ) ;
					++where ;
				}
			}
			return result ;
		}
	}

	CrossDataSetConcordanceComputation::CrossDataSetConcordanceComputation(
		genfile::CohortIndividualSource const& samples,
		std::string const& sample_id_column
	):
		m_samples( samples ),
		m_call_threshhold( 0.9 )
	{
		m_samples.get_column_values( sample_id_column, boost::bind( &SampleIdList::push_back, &m_sample_ids, _2 ) ) ;
	}

	void CrossDataSetConcordanceComputation::set_alternate_dataset( genfile::CohortIndividualSource::UniquePtr samples, genfile::SNPDataSource::UniquePtr snps ) {
		m_alt_dataset_samples = samples ;
		m_alt_dataset_snps = snps ;
		m_alt_dataset_samples->get_column_values( "ID_1", boost::bind( &SampleIdList::push_back, &m_alt_dataset_sample_ids, _2 ) ) ;
		m_sample_mapping = impl::create_sample_mapping( m_sample_ids, m_alt_dataset_sample_ids ) ;
	}

	void CrossDataSetConcordanceComputation::operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, SampleSexes const&, genfile::VariantDataReader&, ResultCallback callback ) {
		using genfile::string_utils::to_string ;
		
		genfile::SNPIdentifyingData alt_snp ;
		bool found = false ;
		while(
			!found
			&&
			m_alt_dataset_snps->get_next_snp_with_specified_position(
				genfile::ignore(),
				genfile::set_value( alt_snp.SNPID() ),
				genfile::set_value( alt_snp.rsid() ),
				genfile::set_value( alt_snp.position().chromosome() ),
				genfile::set_value( alt_snp.position().position() ),
				genfile::set_value( alt_snp.first_allele() ),
				genfile::set_value( alt_snp.second_allele() ),
				snp.get_position()
			)
		) {
			if( alt_snp == snp ) {
				found = true ;
				break ;
			}
		}
		if( found ) {
			genfile::vcf::ThreshholdingGenotypeSetter< Eigen::VectorXd > setter( m_alt_dataset_genotypes, m_call_threshhold, -1, 0, 1, 2 ) ;
			m_alt_dataset_snps->read_variant_data()->get( "genotypes", setter ) ;

			SampleMapping::const_iterator i = m_sample_mapping.begin() ;
			SampleMapping::const_iterator const end_i = m_sample_mapping.end() ;

			for( std::size_t count = 0; i != end_i; ++i, ++count ) {
				double const alt_genotype = m_alt_dataset_genotypes( i->second ) ;
				std::string const stub = to_string( i->first + 1 ) + "~" + to_string( i->second + 1 ) + "(" + m_sample_ids[ i->first ].as< std::string >() + ")" ;
				if( alt_genotype != -1 ) {
					if( genotypes( i->first, alt_genotype ) > m_call_threshhold ) {
						callback( stub + ":concordance", 1 ) ;
					} else {
						callback( stub + ":concordance", 0 ) ;
					}
				} else {
					callback( stub + ":concordance", genfile::MissingValue() ) ;
				}
			}
		}
	}

	std::string CrossDataSetConcordanceComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "CrossDataSetConcordanceComputation" ;
	}
}

