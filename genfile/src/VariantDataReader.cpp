
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <cassert>
#include "genfile/VariantDataReader.hpp"
#include "genfile/vcf/get_set.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"

namespace genfile {

	VariantDataReader::~VariantDataReader() {}

//	VariantDataReader& VariantDataReader::get( std::string const& spec, std::vector< std::vector< Entry > >& data ) {
//		vcf::VectorSetter setter( data ) ;
//		get( spec, setter ) ;
//		return *this ;
//	}

//	VariantDataReader& VariantDataReader::get( std::string const& spec, genfile::SingleSNPGenotypeProbabilities& data ) {
//		vcf::GenotypeSetter< genfile::SingleSNPGenotypeProbabilities > setter( data ) ;
//		get( spec, setter ) ;
//		return *this ;
//	}

	VariantDataReader& VariantDataReader::get( std::string const& spec, PerSampleSetter const& setter ) {
		// Casting away const: this is bad.  But we want this usage for temporary setters.
		PerSampleSetter& nonconst_setter = const_cast< PerSampleSetter& >( setter ) ;
		return get( spec, nonconst_setter ) ;
	}

	VariantDataReader& VariantDataReader::get( std::string const& spec, PerVariantSetter& data ) {
		assert( 0 ) ; // This function should not be called.
	}
}
