
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_VCF_FORMAT_STRICT_METADATA_PARSER_HPP
#define GENFILE_VCF_FORMAT_STRICT_METADATA_PARSER_HPP

#include <map>
#include <string>
#include <boost/noncopyable.hpp>
#include "genfile/vcf/MetadataParser.hpp"

namespace genfile {
	namespace vcf {
		// A strict metadata parser that parses all entries according to (my reading of) the spec.
		struct StrictMetadataParser: public MetadataParser {
			StrictMetadataParser( std::string const& filename ) ;
			StrictMetadataParser( std::string const& spec, std::istream& stream ) ;

			Metadata const& get_metadata() const { return m_metadata ; }
			std::size_t get_number_of_lines() const { return m_metadata.size() ; }

		private:
			std::string const m_spec ;
			std::string m_version ;
			Metadata m_metadata ;
		private:
			void setup( std::string const& spec, std::istream& stream ) ;
			std::string read_version( std::istream& in ) const ;
			Metadata read_metadata( std::istream& in, std::string const& version ) const ;
			bool read_metadata_line( std::istream& in, std::size_t line_number, Metadata* result ) const ;
			std::pair< std::string, std::map< std::string, std::string > > parse_meta_line(
				std::size_t line_number,
				std::string const& line
			) const ;
			std::pair< std::string, std::string > parse_key_value_pair( std::size_t const line_number, std::string const& line ) const ;
			bool validate_meta_key(
				std::size_t line_number,
				std::string const& key
			) const ;
			std::map< std::string, std::string > parse_meta_value(
				std::size_t line_number,
				std::string const& key,
				std::string value
			) const ;
			bool validate_meta_value(
				std::string const& key,
				std::map< std::string, std::string > const& meta_value
			) const ;
		} ;
	}
}

#endif
