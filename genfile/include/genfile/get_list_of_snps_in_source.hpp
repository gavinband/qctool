#ifndef GENFILE_GET_LIST_OF_SNPS_IN_SOURCE_HPP
#define GENFILE_GET_LIST_OF_SNPS_IN_SOURCE_HPP

#include <vector>
#include <boost/function.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"

namespace genfile {
	std::vector< SNPIdentifyingData > get_list_of_snps_in_source(
		SNPDataSource& source,
		boost::function< void( std::size_t, std::size_t ) > progress_callback = boost::function< void( std::size_t, std::size_t ) >() ) ;
}

#endif
