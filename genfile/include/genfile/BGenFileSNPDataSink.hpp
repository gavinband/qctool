#ifndef GENFILE_BGENFILESNPDATASINK_HPP
#define GENFILE_BGENFILESNPDATASINK_HPP

#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/bgen.hpp"
#include "genfile/BasicBGenFileSNPDataSink.hpp"

namespace genfile {
	class BGenFileSNPDataSink: public BasicBGenFileSNPDataSink
	{
	public:
		BGenFileSNPDataSink(
			std::string const& filename,
			std::string const& free_data,
			bgen::uint32_t flags
		)
		: 	BasicBGenFileSNPDataSink( filename, free_data, e_NoCompression, flags )
		{
		}

		~BGenFileSNPDataSink() {
			// We are about to close the file.
			// To write the correct header info, we seek back to the start and rewrite the header block
			// The header comes after the offset which is 4 bytes.
			stream_ptr()->seekp(4, std::ios_base::beg ) ;
			if( stream_ptr()->bad() ) {
				throw FormatUnsupportedError() ;
			}

			write_header_data( *stream_ptr() ) ;
		}
	} ;
}

#endif