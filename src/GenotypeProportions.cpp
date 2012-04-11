
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include "floating_point_utils.hpp"
#include "GenotypeProportions.hpp"

GenotypeProportions::GenotypeProportions()
	: m_proportion_of_AA( 0 ),
	  m_proportion_of_AB( 0 ),
	  m_proportion_of_BB( 0 )
{}

GenotypeProportions::GenotypeProportions( double aa, double ab, double bb )
 : m_proportion_of_AA( aa ),
    m_proportion_of_AB( ab ),
    m_proportion_of_BB( bb )
{}

double& GenotypeProportions::operator[]( std::size_t i ) {
	switch( i ) {
		case(0):
			return m_proportion_of_AA ; break ;
		case(1):
			return m_proportion_of_AB ; break ;
		case(2):
			return m_proportion_of_BB ; break ;
		case(3):
			assert(0) ; break ;
	}
	return m_proportion_of_AA ; //suppress warning.
}

double GenotypeProportions::operator[]( std::size_t i ) const {
	switch( i ) {
		case(0):
			return m_proportion_of_AA ; break ;
		case(1):
			return m_proportion_of_AB ; break ;
		case(2):
			return m_proportion_of_BB ; break ;
		case(3):
			assert(0) ; break ;
	}
	return 0 ; //suppress warning.
}

bool GenotypeProportions::operator==( GenotypeProportions const& other ) const {
	return m_proportion_of_AA == other.m_proportion_of_AA
		&& m_proportion_of_AB == other.m_proportion_of_AB
		&& m_proportion_of_BB == other.m_proportion_of_BB ;
}

bool GenotypeProportions::operator!=( GenotypeProportions const& other ) const {
	return m_proportion_of_AA != other.m_proportion_of_AA
		|| m_proportion_of_AB != other.m_proportion_of_AB
		|| m_proportion_of_BB != other.m_proportion_of_BB ;
}

void GenotypeProportions::floor() {
	m_proportion_of_AA = std::floor( m_proportion_of_AA ) ;
	m_proportion_of_AB = std::floor( m_proportion_of_AB ) ;
	m_proportion_of_BB = std::floor( m_proportion_of_BB ) ;
}

void GenotypeProportions::round() {
	m_proportion_of_AA = round_to_nearest_integer( m_proportion_of_AA ) ;
	m_proportion_of_AB = round_to_nearest_integer( m_proportion_of_AB ) ;
	m_proportion_of_BB = round_to_nearest_integer( m_proportion_of_BB ) ;
}

void GenotypeProportions::zero() {
	m_proportion_of_AA = m_proportion_of_AB = m_proportion_of_BB = 0.0 ;
}

std::ostream& operator<<( std::ostream& aStream, GenotypeProportions const& proportions ) {
	return aStream << proportions.AA() << " " << proportions.AB() << " " << proportions.BB() ;
}

GenotypeProportions operator+( GenotypeProportions const& left, GenotypeProportions const& right ) {
    GenotypeProportions result = left ;
    result += right ;    
    return result ;
}

GenotypeProportions operator/( GenotypeProportions const& left, double right ) {
    GenotypeProportions result = left ;
    result /= right ;    
    return result ;
}

std::vector< GenotypeProportions > operator+( std::vector<GenotypeProportions> const& left, std::vector<GenotypeProportions> const& right ) {
    assert( left.size() == right.size() ) ;
    std::vector<GenotypeProportions > result ( left ) ;
    for( std::size_t i = 0; i < result.size(); ++i ) {
        result[i] += right[i] ;
    }

    return result ;
}

std::vector< GenotypeProportions > operator/( std::vector<GenotypeProportions> const& left, double scalar ) {
    assert( scalar != 0.0 ) ;
    std::vector< GenotypeProportions > result( left ) ;
    for( std::size_t i = 0; i < result.size(); ++i ) {
        result[i] /= scalar ;
    }
    return result ;
}

