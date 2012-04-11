
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <limits>
#include "ConcordanceComputation.hpp"

ConcordanceComputation::ConcordanceComputation( appcontext::OptionProcessor const&, appcontext::UIContext& ):
	m_threshhold( 0.9 )
{}


double ConcordanceComputation::operator()(
	std::size_t const sample1,
	std::size_t const sample2,
	std::vector< genfile::SingleSNPGenotypeProbabilities > const& genotypes
) {
	double non_missing_count = 0.0 ;
	double concordant_count = 0.0 ;
	for( std::size_t snp_i = 0; snp_i < genotypes.size(); ++snp_i ) {
		std::size_t
			genotype1 = std::numeric_limits< std::size_t >::max(),
			genotype2 = std::numeric_limits< std::size_t >::max()
		;

		for( std::size_t g = 0; g < 3; ++g ) {
			if( genotypes[snp_i]( sample1, g ) >= m_threshhold ) {
				genotype1 = g ;
				break ;
			}
		}

		for( std::size_t g = 0; g < 3; ++g ) {
			if( genotypes[snp_i]( sample2, g ) >= m_threshhold ) {
				genotype2 = g ;
				break ;
			}
		}
			
		if( genotype1 < 3 && genotype2 < 3 ) {
			++non_missing_count ;
			if( genotype1 == genotype2 ) {
				++concordant_count ;
			}
		}
	}
	return concordant_count / non_missing_count ;
}

std::string ConcordanceComputation::get_summary() const {
	return "Concordance (proportion of SNPs with identical genotypes) between samples.  Only SNPs well-called in both samples are considered." ;
}

std::string PairwiseNonMissingnessComputation::get_summary() const {
	return "Number of SNPs which are well called in both samples." ;
}

PairwiseNonMissingnessComputation::PairwiseNonMissingnessComputation( appcontext::OptionProcessor const&, appcontext::UIContext& ):
	m_threshhold( 0.9 )
{}

double PairwiseNonMissingnessComputation::operator()(
	std::size_t const sample1,
	std::size_t const sample2,
	std::vector< genfile::SingleSNPGenotypeProbabilities > const& genotypes
) {
	double non_missing_count = 0.0 ;
	for( std::size_t snp_i = 0; snp_i < genotypes.size(); ++snp_i ) {
		std::size_t
			genotype1 = std::numeric_limits< std::size_t >::max(),
			genotype2 = std::numeric_limits< std::size_t >::max()
		;

		for( std::size_t g = 0; g < 3; ++g ) {
			if( genotypes[snp_i]( sample1, g ) >= m_threshhold ) {
				genotype1 = g ;
				break ;
			}
		}

		for( std::size_t g = 0; g < 3; ++g ) {
			if( genotypes[snp_i]( sample2, g ) >= m_threshhold ) {
				genotype2 = g ;
				break ;
			}
		}
			
		if( genotype1 < 3 && genotype2 < 3 ) {
			++non_missing_count ;
		}
	}
	return non_missing_count ;
}

