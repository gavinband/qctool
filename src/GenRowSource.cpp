#include <memory>
#include <boost/bind.hpp>
#include "GenRow.hpp"
#include "GenRowSource.hpp"

SNPDataSourceGenRowSource& SNPDataSourceGenRowSource::read( GenRow & row ){
	m_snp_data_source->read_snp(
		boost::bind< void >( &GenRow::set_number_of_samples, &row, _1 ),
		boost::bind< void >( &GenRow::set_SNPID, &row, _1 ),
		boost::bind< void >( &GenRow::set_RSID, &row, _1 ),
		boost::bind< void >( &GenRow::set_SNP_position, &row, _1 ),
		boost::bind< void >( &GenRow::set_allele1, &row, _1 ),
		boost::bind< void >( &GenRow::set_allele2, &row, _1 ),
		boost::bind< void >( &GenRow::set_genotype_probabilities, &row, _1, _2, _3, _4 )
	) ;
	
	return *this ;
}


