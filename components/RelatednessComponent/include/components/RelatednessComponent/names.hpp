
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_RELATEDNESS_COMPONENT_GET_CONCATENDATED_SAMPLE_IDS_HPP
#define QCTOOL_RELATEDNESS_COMPONENT_GET_CONCATENDATED_SAMPLE_IDS_HPP

#include <vector>
#include <string>
#include "genfile/CohortIndividualSource.hpp"

namespace pca {
	std::string get_concatenated_sample_ids( genfile::CohortIndividualSource const* samples, std::size_t i ) ;
	std::string string_and_number( std::string const& s, std::size_t i ) ;
	std::string get_entry( std::vector< std::string > const&, std::size_t i ) ;
}

#endif
