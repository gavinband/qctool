
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "AlleleProportions.hpp"

AlleleProportions::AlleleProportions( double a, double b )
  : m_proportion_of_A( a ),
    m_proportion_of_B( b )
{}

std::ostream& operator<<( std::ostream& aStream, AlleleProportions const& proportions ) {
    return aStream << "A:" << proportions.proportion_of_A()
        << " B:" << proportions.proportion_of_B() ;
}


GenotypeProportions AlleleProportions::genotype_proportions_at_hardy_weinberg() const {
    GenotypeProportions genotype_proportions( (proportion_of_A() * proportion_of_A()),
												2.0 * (proportion_of_A() * proportion_of_B()),
												(proportion_of_B() * proportion_of_B())) ;
												
	genotype_proportions /= genotype_proportions.sum() ;
	return genotype_proportions ;
}


AlleleProportions operator/( AlleleProportions const& left, double right ) {
	AlleleProportions result = left ;
	return result /= right ;
}
