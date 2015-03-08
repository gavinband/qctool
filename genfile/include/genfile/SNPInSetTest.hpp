
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SNPINSET_TEST_HPP
#define GENFILE_SNPINSET_TEST_HPP

#include <set>
#include <string>
#include "SNPIdentifyingData.hpp"
#include "SNPIdentifyingDataTest.hpp"

namespace genfile {
	struct SNPInSetTest: public SNPIdentifyingDataTest
	{
	public:
		SNPInSetTest( std::set< SNPIdentifyingData > const&, SNPIdentifyingData::CompareFields const& comparer ) ;

		bool operator()( SNPIdentifyingData const& data ) const ;

		bool operator()(
			std::string SNPID,
			std::string RSID,
			GenomePosition position,
			std::string first_allele,
			std::string second_allele
		) const ;
		
		std::string display() const ;

	private:
		std::set< SNPIdentifyingData, SNPIdentifyingData::CompareFields > const m_snps ;
	} ;
}

#endif

