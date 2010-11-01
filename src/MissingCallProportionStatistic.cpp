#include <limits>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <numeric>

#include "GenRow.hpp"
#include "MissingCallProportionStatistic.hpp"

MissingCallProportionStatistic::MissingCallProportionStatistic( double threshhold ):
	m_threshhold( threshhold )
{}

double MissingCallProportionStatistic::calculate_value( GenRowStatistics const& statistics ) const {
	if( statistics.number_of_samples() == 0 ) {
		return std::numeric_limits< double >::quiet_NaN() ;
	}

	GenRow::genotype_proportion_const_iterator
		i( statistics.row().begin_genotype_proportions() ),
		end_i( statistics.row().end_genotype_proportions() ) ;

	double missing = 0.0 ;
	for( ; i != end_i; ++i ) {
		if( std::max( std::max( i->AA(), i->AB() ), i->BB() ) < m_threshhold ) {
			++missing ;
		}
	}
	return missing / statistics.number_of_samples() ;
}
