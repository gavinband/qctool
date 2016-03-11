
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <boost/bind.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/get_list_of_snps_in_source.hpp"
#include "genfile/get_set.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	std::vector< VariantIdentifyingData > get_list_of_snps_in_source(
		SNPDataSource& source,
		boost::function< void( std::size_t, boost::optional< std::size_t > ) > progress_callback
	) {
		typedef std::vector< VariantIdentifyingData > Result ;
		Result result ;
		if( source.total_number_of_snps() ) {
			result.reserve( *source.total_number_of_snps() ) ;
		} else {
			result.reserve( 10000 ) ;
		}
		source.list_snps(
			boost::bind( &Result::push_back, &result, _1 ),
			progress_callback
		) ;
		if( source.total_number_of_snps() && result.size() != *source.total_number_of_snps() ) {
			throw genfile::MalformedInputError( source.get_source_spec(), result.size() ) ;
		}
		return result ;
	}
}
