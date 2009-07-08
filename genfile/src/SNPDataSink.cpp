#include <iostream>
#include <string>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include "snp_data_utils.hpp"
#include "SNPDataSink.hpp"
#include "GenFileSNPDataSink.hpp"
#include "BGenFileSNPDataSink.hpp"

namespace genfile {

	std::auto_ptr< SNPDataSink > SNPDataSink::create(
		std::string const& filename,
		std::string const& free_data
	) {
		return SNPDataSink::create( filename, genfile::filename_indicates_file_is_gzipped( filename ), free_data ) ;
	}

	std::auto_ptr< SNPDataSink > SNPDataSink::create(
		std::string const& filename,
		bool file_is_gzipped,
		std::string const& free_data
	) {
		if( genfile::filename_indicates_bgen_format( filename )) {
			if( file_is_gzipped ) {
				return std::auto_ptr< SNPDataSink >( new ZippedBGenFileSNPDataSink( filename, free_data )) ;
			}
			else {
				return std::auto_ptr< SNPDataSink >( new UnzippedBGenFileSNPDataSink( filename, free_data )) ;
			}
		}
		else {
			return std::auto_ptr< SNPDataSink >( new GenFileSNPDataSink( filename, file_is_gzipped )) ;
		}
	}
}
