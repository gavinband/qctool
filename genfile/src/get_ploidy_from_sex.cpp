
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/CohortIndividualSource.hpp"
#include "genfile/Chromosome.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/get_ploidy_from_sex.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	int get_ploidy_from_sex( CohortIndividualSource const& samples, std::string const& sex_column, Chromosome const& chromosome, std::size_t i ) {
		// This implementation is appropriate for human samples only.
		int result = -1 ;
		if( chromosome == "0X" || chromosome == "0Y" ) {
			VariantEntry entry = samples.get_entry( i, sex_column ) ;
			if( !entry.is_missing() ) {
				std::string const sex = string_utils::to_lower( entry.as< std::string >() ) ;
				if( sex == "1" || sex == "m" || sex == "male" ) {
					result = 1 ;
				} else if( sex == "2" || sex == "f" || sex == "female" ) {
					result = ( chromosome == "0X" ) ? 1 : 0 ;
				} else {
					throw genfile::MalformedInputError( samples.get_source_spec(), "Malformed sex value \"" + sex + "\"", i, samples.get_column_spec().find_column( sex_column ) ) ;
				}
			}
		} else {
			// assume diploid off the sex chromosomes.
			result = 2  ;
		}
		return result ;
	}
}

