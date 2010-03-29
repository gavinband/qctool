#include <iostream>
#include <string>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/GenFileSNPDataSink.hpp"
#include "genfile/BGenFileSNPDataSink.hpp"
#include "genfile/SortingBGenFileSNPDataSink.hpp"

namespace genfile {

	std::auto_ptr< SNPDataSink > SNPDataSink::create(
		std::string const& filename,
		std::string const& free_data,
		bool sort
	) {
		return SNPDataSink::create( filename, get_compression_type_indicated_by_filename( filename ), free_data, sort ) ;
	}

	std::auto_ptr< SNPDataSink > SNPDataSink::create(
		std::string const& filename,
		CompressionType compression_type,
		std::string const& free_data,
		bool sort
	) {
		if( filename_indicates_bgen_format( filename )) {
			assert( compression_type == e_NoCompression ) ;
			if( sort ) {
				return std::auto_ptr< SNPDataSink >( new SortingBGenFileSNPDataSink( filename, free_data, bgen::e_CompressedSNPBlocks )) ;
			}
			else {
				return std::auto_ptr< SNPDataSink >( new BGenFileSNPDataSink( filename, free_data, bgen::e_CompressedSNPBlocks )) ;
			}
		}
		else {
			if( sort ) {
				throw FormatUnsupportedError() ;
			}
			return std::auto_ptr< SNPDataSink >( new GenFileSNPDataSink( filename, get_chromosome_indicated_by_filename( filename ), compression_type )) ;
		}
	}
}
