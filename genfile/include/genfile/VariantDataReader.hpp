#ifndef GENFILE_VARIANT_DATA_READER_HPP
#define GENFILE_VARIANT_DATA_READER_HPP

#include <memory>
#include <vector>
#include <string>
#include "boost/noncopyable.hpp"
#include "genfile/VariantEntry.hpp"

namespace genfile {
	class VariantDataReader: public boost::noncopyable
	{
	public:
		typedef std::auto_ptr< VariantDataReader > UniquePtr ;
		typedef genfile::VariantEntry Entry ;
	public:
		typedef boost::function< void ( std::size_t i, std::vector< Entry > const& ) > Setter ;
	public:
		virtual ~VariantDataReader() {} ;

		virtual VariantDataReader& get( std::string const& spec, Setter setter ) = 0 ;
	} ;
}

#endif
