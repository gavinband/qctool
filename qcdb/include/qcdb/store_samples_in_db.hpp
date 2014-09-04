
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNP_OUTPUT_COMPONENT_STORE_SAMPLES_HPP
#define SNP_OUTPUT_COMPONENT_STORE_SAMPLES_HPP

#include <iostream>
#include <string>
#include "genfile/SNPDataSink.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "qcdb/DBOutputter.hpp"

namespace qcdb {
	void store_samples_in_db(
		std::size_t number_of_samples,
		genfile::SNPDataSink::SampleNameGetter getter,
		genfile::CohortIndividualSource const& sample_data,
		qcdb::DBOutputter& outputter,
		std::string const& table_name
	) ;
}

#endif
