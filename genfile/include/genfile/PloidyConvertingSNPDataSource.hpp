
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_HAPLOID_CORRECTING_SNP_DATA_SOURCE_HPP
#define GENFILE_HAPLOID_CORRECTING_SNP_DATA_SOURCE_HPP

#include <string>
#include <vector>
#include <memory>
#include <boost/function.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/VariantDataReader.hpp"

namespace genfile {
	// A SNPDataSource that wraps an existing source but converts
	// diploid calls / genotype probabilities into haploid calls / probabilities
	// for a specified list of samples.
	// Diploid homozygous calls get converted into haploid calls in the obvious way
	// Diploid heterozygoes calls get converted into missing data.
	struct PloidyConvertingSNPDataSource: public SNPDataSource {
		typedef std::auto_ptr< PloidyConvertingSNPDataSource > UniquePtr ;
		typedef boost::function< int ( genfile::VariantIdentifyingData const&, std::size_t i ) > GetPloidy ;

		static UniquePtr create(
			SNPDataSource::UniquePtr source,
			GetPloidy get_ploidy
		) ;

		PloidyConvertingSNPDataSource(
			SNPDataSource::UniquePtr source,
			GetPloidy get_ploidy
		) ;

	public:
		operator bool() const { return *m_source ; }
		Metadata get_metadata() const { return m_source->get_metadata() ; }
		unsigned int number_of_samples() const { return m_source->number_of_samples() ;}
		void get_sample_ids( GetSampleIds getter ) const { return m_source->get_sample_ids( getter ) ; }

		OptionalSnpCount total_number_of_snps() const { return m_source->total_number_of_snps() ; }
		std::string get_source_spec() const { return "PloidyConverting:" + m_source->get_source_spec() ; }
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
			return prefix + "PloidyConverting:\n" + m_source->get_summary( "  " + prefix, column_width ) ;
		}

	protected:
		void get_snp_identifying_data_impl( VariantIdentifyingData* variant ) {
			m_source->get_snp_identifying_data( &m_variant ) ;
			*variant = m_variant ;
		}
		
		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() {
			m_source->ignore_snp_probability_data() ;
		}

		void reset_to_start_impl() {
			m_source->reset_to_start() ;
		}

	private:
		SNPDataSource::UniquePtr m_source ;
		GetPloidy m_get_ploidy ;
		VariantIdentifyingData m_variant ;
		std::vector< int > m_ploidy ;
	} ;
}

#endif
