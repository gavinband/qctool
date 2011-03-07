#ifndef GENFILE_VARIANTENTRY_HPP
#define GENFILE_VARIANTENTRY_HPP

#include <string>
#include <boost/variant/variant.hpp>
#include "genfile/MissingValue.hpp"

namespace genfile {
	typedef boost::variant< MissingValue, std::string, int, double > EntryData ;
}

#endif
