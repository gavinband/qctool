
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_OVERLAP_HPP
#include <vector>
#include <cassert>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/GenomePosition.hpp"
#include "genfile/overlap.hpp"
#include "genfile/Error.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
	std::vector< std::vector< std::size_t > > get_overlapping_snps(
		std::vector< std::vector< genfile::SNPIdentifyingData > > const& cohort_snps
	)
	{
		std::vector< std::vector< std::size_t > > snp_indices( cohort_snps.size() ) ;
		for( std::size_t i = 0; i < cohort_snps.size(); ++i ) {
			snp_indices[i].resize( cohort_snps[i].size() ) ;
			for( std::size_t j = 0; j < cohort_snps[i].size(); ++j ) {
				snp_indices[i][j] = j ;
			}
		}
		return get_overlapping_snps( cohort_snps, snp_indices ) ;
	}

	std::vector< std::vector< std::size_t > > get_overlapping_snps(
		std::vector< std::vector< genfile::SNPIdentifyingData > > const& cohort_snps,
		std::vector< std::vector< std::size_t > > const& indices_of_snps_to_include
	) {
		assert( indices_of_snps_to_include.size() == cohort_snps.size() ) ;
		
		std::vector< std::vector< std::size_t > > overlapping_snps( cohort_snps.size() ) ;
		std::vector< std::size_t > snp_counter( cohort_snps.size(), 0u ) ;
		
		for( ; snp_counter[0] < indices_of_snps_to_include[0].size(); ++snp_counter[0] ) {
			// In each of the other cohorts, find the next included SNP which has position at least as large
			// as the position of the SNP in cohort 0
			bool reached_end = false ;

			for( std::size_t i = 1; i < cohort_snps.size(); ++i ) {
				// Find the next included SNP in cohort i with position not less than the
				// current SNP in cohort 0.
				for( ;
					snp_counter[i] < indices_of_snps_to_include[i].size() &&
					cohort_snps[i][ indices_of_snps_to_include[i][ snp_counter[i] ] ].get_position()
						< cohort_snps[0][ indices_of_snps_to_include[0][ snp_counter[0] ] ].get_position() ;
					++snp_counter[i]
				) ;

				if( snp_counter[i] == indices_of_snps_to_include[i].size() ) {
					reached_end = true ;
					break ;
				}
			}

			if( reached_end ) {
				break ;
			}

			// Look to see whether this SNP is in the overlap set and is not filtered out
			bool overlap = true ;
			for( std::size_t i = 1; i < cohort_snps.size() && overlap ; ++i ) {
				if( cohort_snps[i][ indices_of_snps_to_include[i][ snp_counter[i] ] ].get_position()
					!= cohort_snps[0][ indices_of_snps_to_include[0][ snp_counter[0] ] ].get_position() ) {
					overlap = false ;
				}
			}
			if( overlap ) {
				for( std::size_t i = 0; i < cohort_snps.size(); ++i ) {
					overlapping_snps[i].push_back( indices_of_snps_to_include[i][ snp_counter[i] ] ) ;
					if( i > 0 ) {
						// make sure we don't look at this SNP again.
						++snp_counter[i] ;
					}
				}
			}
		}
		
		for( std::size_t i = 1; i < cohort_snps.size(); ++i ) {
			assert( overlapping_snps[i].size() == overlapping_snps[0].size() ) ;
		}
		
		return overlapping_snps ;
	}
}

#endif
