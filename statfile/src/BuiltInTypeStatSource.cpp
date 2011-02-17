#include <memory>
#include <string>
#include "statfile/statfile_utils.hpp"
#include "statfile/StatSource.hpp"
#include "statfile/BinFormatStatSource.hpp"
#include "statfile/PackedBinFormatStatSource.hpp"
#include "statfile/RFormatStatSource.hpp"
#include "statfile/DelimitedStatSource.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include "statfile/BuiltInTypeStatSourceChain.hpp"

namespace statfile {
	BuiltInTypeStatSource::UniquePtr BuiltInTypeStatSource::open( std::string const& filename ) {
		FileFormatType format = statfile::get_file_format_type_indicated_by_filename( filename ) ;
		std::auto_ptr< BuiltInTypeStatSource > source ;
		if( format == e_BinFormat ) {
			source.reset( new statfile::BinFormatStatSource( filename )) ;
		}
		else if( format == e_PackedBinFormat ) {
			source.reset( new statfile::PackedBinFormatStatSource( filename )) ;
		}
		else if( format == e_CommaDelimitedFormat ){
			source.reset( new statfile::DelimitedStatSource( filename, "," )) ;
		}
		else if( format == e_TabDelimitedFormat ){
			source.reset( new statfile::DelimitedStatSource( filename, "\t" )) ;
		}
		else if( format == e_RFormat ) {
			source.reset( new statfile::RFormatStatSource( filename )) ;
		}
		else {
			// default to R format
			source.reset( new statfile::RFormatStatSource( filename )) ;
		}

		return source ;
	}
	
	BuiltInTypeStatSource::UniquePtr BuiltInTypeStatSource::open( std::vector< genfile::wildcard::FilenameMatch > const& filenames ) {
		return BuiltInTypeStatSource::UniquePtr(
			BuiltInTypeStatSourceChain::open( filenames ).release()
		) ;
	}
	
	void BuiltInTypeStatSource::read_value( genfile::Chromosome& chromosome ) {
		std::string s ;
		read_value( s ) ;
		chromosome = genfile::Chromosome( s ) ;
	}

	void BuiltInTypeStatSource::read_value( char& c ) {
		std::string s ;
		read_value( s ) ;
		if( s.size() != 1 ) {
			throw genfile::MalformedInputError( "(unknown object of type BuiltInTypeStatSource)", number_of_rows_read() + 1, current_column() + 1 ) ;
		}
		c = s[0] ;
	}
}