#include <iostream>
#include <string>

#include "snp_data_utils.hpp"
#include "SNPDataSource.hpp"
#include "GenFileSNPDataSource.hpp"
#include "BGenFileSNPDataSource.hpp"
#include "SNPDataSourceChain.hpp"

namespace genfile {
	std::auto_ptr< SNPDataSource > SNPDataSource::create( std::string const& filename ) {
		return SNPDataSource::create( filename, genfile::filename_indicates_file_is_gzipped( filename )) ;
	}

	std::auto_ptr< SNPDataSource > SNPDataSource::create( std::string const& filename, bool file_is_gzipped ) {
		if( genfile::filename_indicates_bgen_format( filename )) {
			return std::auto_ptr< SNPDataSource >( new BGenFileSNPDataSource( filename, file_is_gzipped )) ;
		}
		else {
			return std::auto_ptr< SNPDataSource >( new GenFileSNPDataSource( filename, file_is_gzipped )) ;
		}
	}

	std::auto_ptr< SNPDataSource > SNPDataSource::create( std::vector< std::string > const& filenames ) {
		std::auto_ptr< SNPDataSourceChain > chain( new SNPDataSourceChain() ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			chain->add_provider( SNPDataSource::create( filenames[i] )) ;
		}
		return std::auto_ptr< SNPDataSource >( chain.release() ) ;
	}
}
