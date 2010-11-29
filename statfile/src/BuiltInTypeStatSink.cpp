#include <unistd.h>
#include <memory>
#include <string>
#include "statfile/StatSink.hpp"
#include "statfile/RFormatStatSink.hpp"
#include "statfile/BinFormatStatSink.hpp"
#include "statfile/PackedBinFormatStatSink.hpp"

namespace statfile {
	std::auto_ptr< BuiltInTypeStatSink > BuiltInTypeStatSink::open( std::string const& filename ) {
		std::auto_ptr< statfile::BuiltInTypeStatSink > result ;
		FileFormatType format = get_file_format_type_indicated_by_filename( filename ) ;
		
		if( format == statfile::e_BinFormat ) {
			result.reset( new statfile::BinFormatStatSink( filename ) ) ;
		}
		else if( format == statfile::e_PackedBinFormat ) {
			result.reset( new statfile::PackedBinFormatStatSink( filename ) ) ;
		}
		else if( format == statfile::e_RFormat ) {
			result.reset( new statfile::RFormatStatSink( filename ) ) ;
		}
		else {
			// default to R format.
			result.reset( new statfile::RFormatStatSink( filename ) ) ;
		}
		return result ;
	}
	
	std::auto_ptr< BuiltInTypeStatSink > NullBuiltInTypeStatSink::open() {
		return std::auto_ptr< BuiltInTypeStatSink >( new NullBuiltInTypeStatSink ) ;
	}
	
}

