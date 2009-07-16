#include <iostream>
#include <string>

#include "snp_data_utils.hpp"
#include "SNPDataSource.hpp"
#include "GenFileSNPDataSource.hpp"
#include "BGenFileSNPDataSource.hpp"
#include "SNPDataSourceChain.hpp"

namespace genfile {
	std::auto_ptr< SNPDataSource > SNPDataSource::create( std::string const& filename ) {
		return SNPDataSource::create( filename, get_compression_type_indicated_by_filename( filename )) ;
	}

	std::auto_ptr< SNPDataSource > SNPDataSource::create( std::string const& filename, CompressionType compression_type ) {
		if( filename_indicates_bgen_format( filename ) || filename_indicates_bgen_uncompressed_format( filename )) {
			return std::auto_ptr< SNPDataSource >( new BGenFileSNPDataSource( filename, compression_type )) ;
		}
		else {
			return std::auto_ptr< SNPDataSource >( new GenFileSNPDataSource( filename, compression_type )) ;
		}
	}

	std::auto_ptr< SNPDataSource > SNPDataSource::create( std::vector< std::string > const& filenames ) {
		std::auto_ptr< SNPDataSourceChain > chain( new SNPDataSourceChain() ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			chain->add_source( SNPDataSource::create( filenames[i] )) ;
		}
		return std::auto_ptr< SNPDataSource >( chain.release() ) ;
	}
}
