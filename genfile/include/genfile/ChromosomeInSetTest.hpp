
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_CHROMOSOME_IN_SET_TEST_HPP
#define GENFILE_CHROMOSOME_IN_SET_TEST_HPP

#include <set>
#include <string>
#include "genfile/Chromosome.hpp"
#include "genfile/VariantIdentifyingDataTest.hpp"

namespace genfile {
	struct ChromosomeInSetTest: public VariantIdentifyingDataTest
	{
		ChromosomeInSetTest( std::set< Chromosome > const& chromosomes ) ;
		bool operator()( VariantIdentifyingData const& data ) const ;
		std::string display() const ;
	private:
		std::set< Chromosome > m_chromosomes ;
	} ;
}

#endif
