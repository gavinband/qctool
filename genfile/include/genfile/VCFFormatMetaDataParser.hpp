#ifndef GENFILE_VCF_FORMAT_METADATA_PARSER_HPP
#define GENFILE_VCF_FORMAT_METADATA_PARSER_HPP

#include <map>
#include <string>

namespace genfile {
	struct VCFFormatMetaDataParser {
		VCFFormatMetaDataParser( std::string const& spec, std::istream& stream ) ;

		std::map< std::string, std::map< std::string, std::string > > const& get_metadata() const { return m_metadata ;}

	private:
		std::string const m_spec ;
		std::map< std::string, std::map< std::string, std::string > > m_metadata ;
	private:
		std::map< std::string, std::map< std::string, std::string > > read_metadata( std::istream& in ) const ;
		std::pair< std::string, std::map< std::string, std::string > > parse_meta_line(
			std::size_t line_number,
			std::string const& line
		) const ;
		std::pair< std::string, std::string > parse_key_value_pair( std::string const& line ) const ;
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
		VCFFormatMetaDataParser( VCFFormatMetaDataParser const& other ) ;
	} ;
}

#endif
