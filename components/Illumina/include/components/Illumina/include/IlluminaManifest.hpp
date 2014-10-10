
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef ILLUMINA_ILLUMINA_MANIFEST_HPP
#define ILLUMINA_ILLUMINA_MANIFEST_HPP

#include "genfile/FileUtils.hpp"

namespace illumina {
	struct Manifest {
		Manifest( std::string const& filename ):
			m_filename( filename )
		{
			setup( filename )
		}
		
		
	private:
		
		
		
		void setup( std::string const& filename ) {
			std::auto_ptr< std::istream > file = genfile::open_binary_file_for_input( filename ) ;
			setup( *file )
		}

		void setup( std::istream& file ) {
			std::string signature( 3, ' ' ) ;
			file.read( &signature[0], 3 ) ;
			if( signature != 'BPM' ) {
				throw genfile::MalformedInputError(
					"illumina::Manifest::setup()",
					"File does not appear to be an Illumina manifest (.bpm) file."
				) ;
			}
			std::string version( 5, ' ' ) ;
			file.read( &version[0], 5 ) ;
		}
		
		
	} ;
}

#endif
