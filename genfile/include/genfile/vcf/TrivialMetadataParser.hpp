
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_VCF_TRIVIAL_METADATA_PARSER_HPP
#define GENFILE_VCF_TRIVIAL_METADATA_PARSER_HPP


#include "genfile/vcf/MetadataParser.hpp"

namespace genfile {
	namespace vcf {
		// A metadata parser that ignores all the metadata (lines beginning with ##) and returns an empty map.
		struct TrivialMetadataParser: public MetadataParser {
			TrivialMetadataParser( std::string const& filename ) ;
			TrivialMetadataParser( std::string const& spec, std::istream& stream ) ;

			Metadata const& get_metadata() const { return m_metadata ; }
			std::size_t get_number_of_lines() const { return m_number_of_lines ; }
		private:
			std::string const m_spec ;
			std::size_t m_number_of_lines ;
			Metadata const m_metadata ;
			// Read through lines until either a line not starting with ##, or EOF is found.
			void setup( std::istream& stream ) ;
			bool read_metadata_line( std::istream& in, std::size_t line_number ) const ;
		} ;
	}
}

#endif
