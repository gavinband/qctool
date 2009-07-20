#include "GenotypeProportions.hpp"
#include "ExternalStorageGenRow.hpp"

ExternalStorageGenRow::ExternalStorageGenRow( GenotypeProportions* storage, std::size_t storage_size )
	: m_number_of_samples( storage_size )
{
	assert( storage ) ;
	set_storage( storage, storage_size ) ;
}

void ExternalStorageGenRow::set_storage( GenotypeProportions* storage, std::size_t storage_size )
{
	assert( storage ) ;
	assert( m_number_of_samples <= storage_size ) ;
	m_storage = storage ;
	m_storage_size = storage_size ;
}

void ExternalStorageGenRow::set_number_of_samples( std::size_t n ) {
	assert( n <= m_storage_size ) ;
	m_number_of_samples = n ;
}

void ExternalStorageGenRow::add_genotype_proportions( GenotypeProportions const& genotype_proportions ) {
	assert( m_number_of_samples < m_storage_size ) ;
	*(m_storage + m_number_of_samples) = genotype_proportions ;
}

