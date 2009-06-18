#ifndef __GTOOL_ALLELEPROPORTIONS_HPP__
#define __GTOOL_ALLELEPROPORTIONS_HPP__

#include <numeric>
#include <iostream>
#include "GenotypeProportions.hpp"

struct AlleleProportions {

		// Constructors and assignment
        AlleleProportions( double a, double b ) ;
		AlleleProportions( AlleleProportions const& other )
			: m_proportion_of_A( other.m_proportion_of_A ),
				m_proportion_of_B( other.m_proportion_of_B )
		{}

		AlleleProportions& operator=( AlleleProportions const& other ) {
			m_proportion_of_A = other.m_proportion_of_A ;
			m_proportion_of_B = other.m_proportion_of_B ;
		}
		

		// data access methods
        double proportion_of_A() const { return m_proportion_of_A ; }
        double proportion_of_B() const { return m_proportion_of_B ; }
        double A() const { return m_proportion_of_A ; }
        double B() const { return m_proportion_of_B ; }
        double minor() const { return std::min( m_proportion_of_A, m_proportion_of_B ) ; }
        double major() const { return std::max( m_proportion_of_A, m_proportion_of_B ) ; }
        double sum() const { return m_proportion_of_A + m_proportion_of_B ; }

        GenotypeProportions genotype_proportions_at_hardy_weinberg() const ;

		
		// useful operators
		AlleleProportions& operator/=( double other ) { m_proportion_of_A /= other ; m_proportion_of_B /= other ; return *this ; }

	private:

        double m_proportion_of_A ;
        double m_proportion_of_B ;
} ;

// useful operator
AlleleProportions operator/( AlleleProportions const& left, double right ) ;

std::ostream& operator<<( std::ostream&, AlleleProportions const& ) ;

#endif

