#include "GenRow.hpp"
#include "fill_genotype_probabilities.hpp"

void fill_genotype_probabilities( GenRow* row, double fill_AA, double fill_AB, double fill_BB ) {
	GenRow::genotype_proportion_iterator i = row->begin_genotype_proportions(),
		end_i = row->end_genotype_proportions() ;
	for( ; i != end_i; ++i ) {
		double sum = i->sum() ;
		assert( sum >= 0.0 ) ;
		i->AA() += (1-sum) * fill_AA ;
		i->AB() += (1-sum) * fill_AB ;
		i->BB() += (1-sum) * fill_BB ;
	}
}

