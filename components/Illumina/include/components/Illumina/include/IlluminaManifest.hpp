
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
			setup( *file, filename )
		}

		void setup( std::istream& file, std::string const& filename ) {
			std::string signature( 3, ' ' ) ;
			file.read( &signature[0], 3 ) ;
			if( signature != 'BPM' ) {
				throw genfile::MalformedInputError(
					"illumina::Manifest::setup()",
					"File (" + filename + ") does not appear to be an Illumina manifest (.bpm) file."
				) ;
			}
			unsigned char version_major ;
			uint32_t version_minor ;
			file.read( &version_major, 1 ) ;
			file.read( &version_minor, 4 ) ;
			if( version_major != 1 || version_minor != 4 ) {
				throw genfile::MalformedInputError(
					"illumina::Manifest::setup()",
					"Illumina manifest (" + filename + ") does not appear to have the right version (" + int( version_major ) + "," + version_minor + " instead of 1,4)"
				) ;
			}

			
		}
	} ;
}

#endif
