#include <vector>
#include <boost/bind.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/get_list_of_snps_in_source.hpp"
#include "genfile/get_set.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	std::vector< SNPIdentifyingData > get_list_of_snps_in_source(
		SNPDataSource& source,
		boost::function< void( std::size_t, std::size_t ) > progress_callback
	) {
		typedef std::vector< SNPIdentifyingData > Result ;
		Result result ;
		result.reserve( source.total_number_of_snps() ) ;
		source.list_snps(
			boost::bind( &Result::push_back, &result, _1 ),
			progress_callback
		) ;
		if( result.size() != source.total_number_of_snps() ) {
			throw genfile::MalformedInputError( source.get_source_spec(), result.size() ) ;
		}
		return result ;
	}
}
