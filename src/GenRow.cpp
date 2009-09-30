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
#include "floating_point_utils.hpp"

bool GenRowIdentifyingData::operator==( GenRowIdentifyingData const& right ) const {
	return (
		( m_SNPID == right.m_SNPID )
		&& ( m_RSID == right.m_RSID )
		&& ( m_chromosome == right.m_chromosome )
		&& ( m_SNP_position == right.m_SNP_position )
		&& ( m_1st_allele == right.m_1st_allele )
		&& ( m_2nd_allele == right.m_2nd_allele )
	) ;
}

bool GenRow::operator==( GenRow const& right ) const {
	if(
		GenRowIdentifyingData::operator==( right )
		&& ( number_of_samples() == right.number_of_samples() )
	) {
		genotype_proportion_const_iterator
			i( begin_genotype_proportions() ),
			end_i( end_genotype_proportions() ),
			j( right.begin_genotype_proportions() ) ;
		for( ; i != end_i; ++i, ++j ) {
			if( *i != *j ) {
				return false ;
			}
		} 
		return true ;
	}

	return false ;
}

bool GenRow::operator!=( GenRow const& right ) const {
	return !(*this == right) ;
}

std::size_t GenRow::number_of_samples() const {
	return std::distance( begin_genotype_proportions(), end_genotype_proportions() ) ; 
}

bool GenRow::check_if_equal( GenRow const& left, GenRow const& right, double epsilon ) {
	if( ! (static_cast< GenRowIdentifyingData >( left ) == right )) {
		return false ;
	}
	if( left.number_of_samples() != right.number_of_samples() ) {
		return false ;
	}
	
	genotype_proportion_const_iterator
		i = left.begin_genotype_proportions(),
		end_i = left.end_genotype_proportions(),
		j = right.begin_genotype_proportions() ;
		
	for( ; i != end_i ; ++i, ++j ) {
		if( !floats_are_equal_to_within_epsilon( i->AA(), j->AA(), epsilon )) {
			return false ;
		}
		else if ( !floats_are_equal_to_within_epsilon( i->AB(), j->AB(), epsilon )) {
			return false ;
		}
		else if( !floats_are_equal_to_within_epsilon( i->BB(), j->BB(), epsilon )) {
			return false ;
		}
	}
	
	return true ;
}

void InternalStorageGenRow::filter_out_samples_with_indices( std::vector< std::size_t > const& indices_to_filter_out ) {
	if( indices_to_filter_out.empty()) {
		return ;
	}
	assert( indices_to_filter_out.size() <= m_genotype_proportions.size() ) ;

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
		assert( std::distance( g_1st, g_2nd ) >= 0 ) ;
		std::copy( g_1st, g_2nd, std::back_inserter( new_genotype_proportions )) ;
	}

	g_2nd = m_genotype_proportions.end() ;
	std::copy( g_1st, g_2nd, std::back_inserter( new_genotype_proportions )) ;

	assert( new_genotype_proportions.size() == m_genotype_proportions.size() - indices_to_filter_out.size() ) ;

	m_genotype_proportions = new_genotype_proportions ;
}

