#ifndef __GTOOL_GENROWFILESOURCE_HPP
#define __GTOOL_GENROWFILESOURCE_HPP

#include <string>
#include <vector>
#include <stddef.h>
#include <boost/bind.hpp>

#include "GenRow.hpp"
#include "SimpleFileObjectSource.hpp"
#include "ChainingFileObjectSource.hpp"
#include "FileUtil.hpp"
#include "SNPDataSource.hpp"


// This class is an adapter between SNPDataSource and ObjectSource< GenRow >.
struct SNPDataSourceGenRowSource: public ObjectSource< GenRow >
{
	SNPDataSourceGenRowSource( std::auto_ptr< genfile::SNPDataSource > snp_data_source ) 
		: m_snp_data_source( snp_data_source )
	{}

	SNPDataSourceGenRowSource& read( GenRow & row ) ;

	bool fail() const { return false ; } // TODO: implement this.
	operator bool() { return *m_snp_data_source ; }
	unsigned int number_of_samples() const { return m_snp_data_source->number_of_samples() ; }
	unsigned int total_number_of_snps() const { return m_snp_data_source->total_number_of_snps() ; }

private:

	std::auto_ptr< genfile::SNPDataSource > m_snp_data_source ;
} ;

#endif
