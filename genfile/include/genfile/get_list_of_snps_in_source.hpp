
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_GET_LIST_OF_SNPS_IN_SOURCE_HPP
#define GENFILE_GET_LIST_OF_SNPS_IN_SOURCE_HPP

#include <vector>
#include <boost/function.hpp>
#include <boost/optional.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/VariantIdentifyingData.hpp"

namespace genfile {
	std::vector< VariantIdentifyingData > get_list_of_snps_in_source(
		SNPDataSource& source,
		boost::function< void( std::size_t, boost::optional< std::size_t > ) > progress_callback
			= boost::function< void( std::size_t, boost::optional< std::size_t > ) >()
	) ;
}

#endif
