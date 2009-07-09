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
		return SNPDataSink::create( filename, get_compression_type_indicated_by_filename( filename ), free_data ) ;
	}

	std::auto_ptr< SNPDataSink > SNPDataSink::create(
		std::string const& filename,
		CompressionType compression_type,
		std::string const& free_data
	) {
		if( filename_indicates_bgen_format( filename )) {
			if( compression_type == e_GzipCompression ) {
				return std::auto_ptr< SNPDataSink >( new ZippedBGenFileSNPDataSink( filename, free_data )) ;
			}
			else {
				return std::auto_ptr< SNPDataSink >( new BGenFileSNPDataSink( filename, free_data, bgen::e_NoFlags )) ;
			}
		}
		else if( filename_indicates_bgen_compressed_format( filename )) {
			return std::auto_ptr< SNPDataSink >( new BGenFileSNPDataSink( filename, free_data, bgen::e_CompressedSNPBlocks )) ;			
		}
		else {
			return std::auto_ptr< SNPDataSink >( new GenFileSNPDataSink( filename, compression_type )) ;
		}
	}
}
