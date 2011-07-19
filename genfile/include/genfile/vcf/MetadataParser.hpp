#ifndef GENFILE_VCF_FORMAT_METADATA_PARSER_HPP
#define GENFILE_VCF_FORMAT_METADATA_PARSER_HPP

#include <map>
#include <string>
#include <memory>
#include <boost/noncopyable.hpp>

namespace genfile {
	namespace vcf {
		struct MetadataParser: public boost::noncopyable {
			typedef std::auto_ptr< MetadataParser > UniquePtr ;
			virtual ~MetadataParser() {} ;
			typedef std::multimap< std::string, std::map< std::string, std::string > > Metadata ;
			virtual Metadata const& get_metadata() const = 0 ;
			virtual std::size_t get_number_of_lines() const = 0 ;
		} ;
	}
}

#endif
