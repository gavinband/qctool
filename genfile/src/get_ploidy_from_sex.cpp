
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
	int get_ploidy_from_sex(
		std::vector< char > const& sexes,
		genfile::VariantIdentifyingData const& variant,
		std::size_t sample_i
	) {
		genfile::Chromosome const& chromosome = variant.get_position().chromosome() ;
		if( chromosome == "0X" || chromosome == "X" || chromosome == "chrX" ) {
			if( sexes[sample_i] == 'm' ) {
				return 1 ;
			} else if( sexes[sample_i] == 'f' ) {
				return 2 ; 
			} else {
				return -1 ;
			}
		} else if( chromosome == "0Y" || chromosome == "Y" || chromosome == "chrY" ) {
			if( sexes[sample_i] == 'm' ) {
				return 1 ;
			} else if( sexes[sample_i] == 'f' ) {
				return 0 ; 
			} else {
				return -1 ;
			}
		} else {
			return 2 ;
		}
	}
}

