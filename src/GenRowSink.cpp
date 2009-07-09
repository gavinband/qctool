#include <memory>
#include <boost/bind.hpp>
#include "GenRowSink.hpp"
#include "SNPDataSink.hpp"

SNPDataSinkGenRowSink::SNPDataSinkGenRowSink( std::auto_ptr< genfile::SNPDataSink > snp_data_sink )
: m_snp_data_sink( snp_data_sink )
{}

SNPDataSinkGenRowSink& SNPDataSinkGenRowSink::write( GenRow const& row ) {
	m_snp_data_sink->write_snp(
		row.number_of_samples(),
		row.SNPID(),
		row.RSID(),
		row.SNP_position(),
		row.first_allele(),
		row.second_allele(),
		boost::bind< double >( &GenRow::get_AA_probability, &row, _1 ),
		boost::bind< double >( &GenRow::get_AB_probability, &row, _1 ),
		boost::bind< double >( &GenRow::get_BB_probability, &row, _1 )
	) ;
	
	return *this ;
}

