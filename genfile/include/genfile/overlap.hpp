
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_OVERLAP_HPP
#include <vector>
#include "genfile/VariantIdentifyingData.hpp"

namespace genfile {

	// Given a list of SNPs in each cohort, return a vector of sets of indices, one set per cohort.
	// The ith set consists of indices of SNPs in cohort i that are in the overlap set with SNPs in all other cohorts.
	// Precondition: each of the lists in cohort_snps must store SNPs in nondecreasing position order.
	// Postcondition: all result sets have the same size.
	std::vector< std::vector< std::size_t > > get_overlapping_snps(
		std::vector< std::vector< genfile::VariantIdentifyingData > > const& cohort_snps
	) ;
	

	// As above but ignore any SNP in cohort i whose index is not listed in the corresponding entry of the second argument.
	// (The entries of the second argument are assumed to be in nondecreasing order).
	std::vector< std::vector< std::size_t > > get_overlapping_snps(
		std::vector< std::vector< genfile::VariantIdentifyingData > > const& cohort_snps,
		std::vector< std::vector< std::size_t > > const& indices_of_snps_to_include
	) ;
}

#endif
