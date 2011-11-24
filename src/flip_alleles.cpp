#ifndef GENFILE_FLIP_ALLELES_HPP
#define GENFILE_FLIP_ALLELES_HPP

#include "GenRow.hpp"
#include "flip_alleles.hpp"

void flip_alleles( GenRow* row ) {
	std::string const old_second_allele = row->second_allele() ;
	row->set_allele2( row->first_allele() ) ;
	row->set_allele1( old_second_allele ) ;

	GenRow::genotype_proportion_iterator
		i = row->begin_genotype_proportions(),
		end_i = row->end_genotype_proportions() ;
	for( ; i != end_i; ++i ) {
		double BB = i->BB() ;
		i->BB() = i->AA() ;
		i->AA() = BB ;
	}
}

#endif
