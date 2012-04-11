
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GTOOL_EXTERNALSTORAGEGENROW_HPP__
#define GTOOL_EXTERNALSTORAGEGENROW_HPP__

#include "GenotypeProportions.hpp"
#include "GenRow.hpp"

// A GenRow which holds its GenotypeProportions in externally-specified
// storage.
// Invariant: number_of_samples() <= size of storage
class ExternalStorageGenRow: public GenRow
{
public:

	ExternalStorageGenRow( GenotypeProportions* storage, std::size_t storage_size ) ;
	void set_storage( GenotypeProportions* storage, std::size_t storage_size ) ;

	void set_number_of_samples( std::size_t n ) ;

	genotype_proportion_const_iterator begin_genotype_proportions() const { return m_storage ; }
	genotype_proportion_const_iterator end_genotype_proportions() const { return m_storage + m_number_of_samples ; }
	genotype_proportion_iterator begin_genotype_proportions() { return m_storage ; }
	genotype_proportion_iterator end_genotype_proportions() { return m_storage + m_number_of_samples ; }

	void add_genotype_proportions( GenotypeProportions const& genotype_proportions ) ;

private:
	GenotypeProportions* m_storage ;
	std::size_t m_storage_size ;
	std::size_t m_number_of_samples ;
} ;




#endif

