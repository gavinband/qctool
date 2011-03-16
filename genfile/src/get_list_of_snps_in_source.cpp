#include <vector>
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
		std::vector< SNPIdentifyingData > result ;
		result.reserve( source.total_number_of_snps() ) ;
		source.reset_to_start() ;
		for(
			SNPIdentifyingData data ;
			source.get_snp_identifying_data(
				genfile::ignore(),
				genfile::set_value( data.SNPID() ),
				genfile::set_value( data.rsid() ),
				genfile::set_value( data.position().chromosome() ),
				genfile::set_value( data.position().position() ),
				genfile::set_value( data.first_allele() ),
				genfile::set_value( data.second_allele() )
			) ;
			source.ignore_snp_probability_data()
		) {
			result.push_back( data ) ;
			if( progress_callback ) {
				progress_callback( source.number_of_snps_read() + 1, source.total_number_of_snps() ) ;
			}
		}
		if( result.size() != source.total_number_of_snps() ) {
			throw genfile::MalformedInputError( source.get_source_spec(), result.size() ) ;
		}
		return result ;
	}
}
