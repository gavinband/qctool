#ifndef QCTOOL_RELATOTRON_HPP
#define QCTOOL_RELATOTRON_HPP

#include <limits>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"

#include "OstreamTee.hpp"

struct Relatotron: public genfile::SNPDataSourceProcessor::Callback
{
	typedef genfile::SingleSNPGenotypeProbabilities SingleSNPGenotypeProbabilities ;
	typedef genfile::SNPIdentifyingData SNPIdentifyingData ;

	Relatotron( OstreamTee& logger ):
	 	m_logger( &logger )
	{} ;
	
	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
		m_number_of_samples = number_of_samples ;
		m_snps.clear() ;
		m_genotypes.clear() ;
		m_allele_frequencies.clear() ;
		m_snps.reserve( number_of_snps ) ;
		m_genotypes.reserve( number_of_snps ) ;
		m_allele_frequencies.reserve( number_of_snps ) ;		
	}

	void processed_snp( SNPIdentifyingData const& id_data, SingleSNPGenotypeProbabilities const& genotypes ) {
		assert( genotypes.get_number_of_samples() == m_number_of_samples ) ;
		assert( m_snps.size() == m_genotypes.size() ) ;
		m_snps.push_back( id_data ) ;
		m_genotypes.push_back( genotypes ) ;
		m_allele_frequencies.push_back( estimate_allele_frequency( genotypes )) ;
	}

	void end_processing_snps() {
		assert( m_snps.size() == m_genotypes.size() ) ;
		assert( m_snps.size() == m_allele_frequencies.size() ) ;

		std::size_t estimated_memory_usage = ( m_snps.capacity() * sizeof( SNPIdentifyingData ) )
			+ ( m_genotypes.capacity() * sizeof( SingleSNPGenotypeProbabilities ))
			+ ( m_genotypes.capacity() * m_number_of_samples * 3 * sizeof( double ))
			+ m_allele_frequencies.size() * sizeof( double ) ;
		
		(*m_logger)
			<< "Relatotron: finished loading "
			<< m_snps.size()
			<< " SNPs, memory usage is "
			<< std::fixed << std::setprecision( 1 ) << ( estimated_memory_usage / 1000000.0 )
			<< "Mb.\n" ;
			
		compute_pairwise_relatedness_coefficients() ;
	}
	
private:
	
	OstreamTee* m_logger ;
	std::size_t m_number_of_samples ;
	
	std::vector< SNPIdentifyingData > m_snps ;
	std::vector< SingleSNPGenotypeProbabilities > m_genotypes ;
	std::vector< double > m_allele_frequencies ;
	
private:
	double estimate_allele_frequency( SingleSNPGenotypeProbabilities const& genotypes ) const {
		// We only use the calls with sum > 0.95 to estimate the allele frequency.
		// I scale those genotypes to have sum 1.  All others are ignored.
		// If there are less than 100 genotypes to estimate with, I set the allele frequency to NaN.
		//
		// Ideally we should use prior information to fill in the missing genotypes.
		// However, that involves knowing the frequencies in e.g. HapMap
		// and the strand alignment.
		double allele_count = 0.0 ;
		std::size_t sample_count = 0 ;

		for( std::size_t i = 0; i < genotypes.size(); ++i ) {
			if( genotypes.sum( i ) >= 0.95 ) {
				allele_count += ( genotypes.AB( i ) + 2.0 * genotypes.BB( i ) ) / genotypes.sum( i ) ;
				++sample_count ;
			}
		}
		
		if( sample_count < 100 ) {
			return std::numeric_limits< double >::quiet_NaN() ;
		}
		else {
			return allele_count / sample_count ;
		}
	}

	void compute_pairwise_relatedness_coefficients() const {
		for( std::size_t sample1 = 0; sample1 < m_number_of_samples; ++sample1 ) {
			(*m_logger) << "Running sample " << (sample1+1) << " of " << (m_number_of_samples ) << "...\n" ;
			for( std::size_t sample2 = sample1; sample2 < m_number_of_samples; ++sample2 ) {
				for( std::size_t snp_i = 0; snp_i < m_snps.size(); ++snp_i ) {
					// do relatedness here.
				}
			}
		}
	}
} ;

#endif
