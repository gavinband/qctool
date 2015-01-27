
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

	MissingnessHeterozygosityComputation::MissingnessHeterozygosityComputation( genfile::Chromosome chromosome ):
		m_chromosome( chromosome ),
	 	m_snp_index( 0 ),
		m_threshhold( 0.9 )
	{}

	void MissingnessHeterozygosityComputation::accumulate( genfile::SNPIdentifyingData const& snp, Genotypes const& genotypes, genfile::VariantDataReader& ) {
		if( m_chromosome && m_chromosome.get() != snp.get_position().chromosome() ) {
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

	void MissingnessHeterozygosityComputation::compute( int sample, ResultCallback callback ) {
		callback( sample, "missing_proportion", ( m_snp_index - m_total_probabilities( sample ) ) / m_snp_index ) ;
		callback( sample, "missing_call_proportion", ( m_snp_index - m_total_calls( sample ) ) / m_snp_index ) ;
		callback( sample, "heterozygous_proportion", m_het_snps( sample ) / m_total_probabilities( sample ) ) ;
		callback( sample, "heterozygous_call_proportion", m_het_snp_calls( sample ) / m_total_calls( sample ) ) ;
	}

	std::string MissingnessHeterozygosityComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "MissingnessHeterozygosityComputation" ;
	}
}

