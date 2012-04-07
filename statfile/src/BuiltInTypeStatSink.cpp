#include <unistd.h>
#include <memory>
#include <string>
#include "statfile/StatSink.hpp"
#include "statfile/RFormatStatSink.hpp"
#include "statfile/BinFormatStatSink.hpp"
#include "statfile/PackedBinFormatStatSink.hpp"
#include "statfile/DelimitedStatSink.hpp"

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
		else if( format == statfile::e_CommaDelimitedFormat ) {
			result.reset( new statfile::DelimitedStatSink( filename, "," ) ) ;
		}
		else if( format == statfile::e_TabDelimitedFormat ) {
			result.reset( new statfile::DelimitedStatSink( filename, "\t" ) ) ;
		}
		else if( format == statfile::e_RFormat ) {
			result.reset( new statfile::DelimitedStatSink( filename, " " ) ) ;
		}
		else {
			// default to tab-separated format.
			result.reset( new statfile::DelimitedStatSink( filename, "\t" ) ) ;
		}
		return result ;
	}
	
	std::auto_ptr< BuiltInTypeStatSink > NullBuiltInTypeStatSink::open() {
		return std::auto_ptr< BuiltInTypeStatSink >( new NullBuiltInTypeStatSink ) ;
	}
}

