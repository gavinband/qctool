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

bool GenRow::operator==( GenRow const& right ) const {
	return ( m_SNPID == right.m_SNPID )
		&& ( m_RSID == right.m_RSID )
		&& ( m_SNP_position == right.m_SNP_position )
		&& ( m_1st_allele == right.m_1st_allele )
		&& ( m_2nd_allele == right.m_2nd_allele )
		&& ( m_genotype_proportions == right.m_genotype_proportions ) ;

}

void GenRow::filter_out_samples_with_indices( std::vector< std::size_t > const& indices_to_filter_out ) {
	if( indices_to_filter_out.empty()) {
		return ;
	}
	assert( indices_to_filter_out.size() < m_genotype_proportions.size() ) ;

	std::vector< GenotypeProportions > new_genotype_proportions ;
	new_genotype_proportions.reserve( m_genotype_proportions.size() - indices_to_filter_out.size() ) ;

	std::vector< std::size_t >::const_iterator i = indices_to_filter_out.begin() ;
	assert( *i < m_genotype_proportions.size() ) ;
	std::vector< GenotypeProportions >::const_iterator
		g_1st = m_genotype_proportions.begin(),
	 	g_2nd = m_genotype_proportions.begin() + *i ;
	std::copy( g_1st, g_2nd, std::back_inserter( new_genotype_proportions )) ;
	
	for( ++i, g_1st = g_2nd + 1; i != indices_to_filter_out.end(); ++i, g_1st = g_2nd + 1 ) {
		assert( *i < m_genotype_proportions.size() ) ;
		g_2nd = m_genotype_proportions.begin() + *i ;
		assert( std::distance( g_1st, g_2nd ) > 0 ) ;
		std::copy( g_1st, g_2nd, std::back_inserter( new_genotype_proportions )) ;
	}

	g_2nd = m_genotype_proportions.end() ;
	std::copy( g_1st, g_2nd, std::back_inserter( new_genotype_proportions )) ;

	assert( new_genotype_proportions.size() == m_genotype_proportions.size() - indices_to_filter_out.size() ) ;

	m_genotype_proportions = new_genotype_proportions ;
}
