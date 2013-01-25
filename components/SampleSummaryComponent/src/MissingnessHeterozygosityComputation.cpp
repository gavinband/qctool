
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <Eigen/Core>
#include "genfile/vcf/get_set_eigen.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"
#include "components/SampleSummaryComponent/MissingnessHeterozygosityComputation.hpp"

namespace sample_stats {

	MissingnessHeterozygosityComputation::MissingnessHeterozygosityComputation():
	 	m_snp_index( 0 ),
		m_threshhold( 0.9 )
	{}

	void MissingnessHeterozygosityComputation::accumulate( genfile::SNPIdentifyingData const& snp, Genotypes const& genotypes, genfile::VariantDataReader& ) {
		// ignore the X chromosome for now.
		if( snp.get_position().chromosome().is_sex_determining() ) {
			return ;
		}

		if( m_snp_index == 0 ) {
			m_total_probabilities.setZero( genotypes.rows() ) ;
			m_total_calls.setZero( genotypes.rows() ) ;
			m_het_snps.setZero( genotypes.rows() ) ;
			m_het_snp_calls.setZero( genotypes.rows() ) ;
			m_ones.setOnes( genotypes.rows() ) ;
		}

		m_total_probabilities += genotypes.rowwise().sum() ;
		m_total_calls.array() += ( genotypes.rowwise().maxCoeff().array() >= m_threshhold ).cast< double >() ;

		m_het_snps += genotypes.col( 1 ) ;
		m_het_snp_calls.array() += ( genotypes.col( 1 ).array() >= m_threshhold ).cast< double >() ;

		m_snp_index++ ;
	}

	void MissingnessHeterozygosityComputation::compute( ResultCallback callback ) {
		for( int sample = 0; sample < m_total_probabilities.size(); ++sample ) {
			callback( sample, "missing proportion", ( m_snp_index - m_total_probabilities( sample ) ) / m_snp_index ) ;
			callback( sample, "missing call proportion", ( m_snp_index - m_total_calls( sample ) ) / m_snp_index ) ;
			callback( sample, "heterozygous proportion", m_het_snps( sample ) / m_total_probabilities( sample ) ) ;
			callback( sample, "heterozygous call proportion", m_het_snp_calls( sample ) / m_total_calls( sample ) ) ;
		}
	}

	std::string MissingnessHeterozygosityComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "MissingnessHeterozygosityComputation" ;
	}
}

