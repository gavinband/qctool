#include <sstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cmath>

#include "../config.hpp"
#include "GToolException.hpp"
#include "GenRow.hpp"
#include "Whitespace.hpp"
#include "GenotypeProportions.hpp"

std::size_t GenRow::number_of_columns() const {
	return (m_genotype_proportions.size() * 3) + std::size_t(5) ;
}


std::size_t GenRow::number_of_samples() const {
	return m_genotype_proportions.size() ;
}

GenotypeProportions const& GenRow::genotype_proportions_for_sample( std::size_t sampleNumber ) const {
	assert( sampleNumber < number_of_samples() ) ;
	return m_genotype_proportions[ sampleNumber ] ;
}

void GenRow::reserveSpaceForNSamples( std::size_t number_of_samples ) {
	m_genotype_proportions.reserve( number_of_samples ) ;
}

bool GenRow::operator==( GenRow const& right ) const {
	return ( m_SNPID == right.m_SNPID )
		&& ( m_RSID == right.m_RSID )
		&& ( m_SNP_position == right.m_SNP_position )
		&& ( m_1st_allele == right.m_1st_allele )
		&& ( m_2nd_allele == right.m_2nd_allele )
		&& ( m_genotype_proportions == right.m_genotype_proportions ) ;

}

