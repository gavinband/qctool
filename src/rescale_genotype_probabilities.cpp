
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "GenRow.hpp"
#include "rescale_genotype_probabilities.hpp"

void rescale_genotype_probabilities( GenRow* row, double ignore_threshhold ) {
	GenRow::genotype_proportion_iterator i = row->begin_genotype_proportions(),
		end_i = row->end_genotype_proportions() ;
	for( ; i != end_i; ++i ) {
		double sum = i->sum() ;
		assert( sum >= 0.0 ) ;
		if( sum < ignore_threshhold ) {
			i->AA() = i->AB() = i->BB() = 0.0 ;
		}
		else {
			(*i) /= sum ;
		}
	}
}
