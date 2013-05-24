
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include <memory>
#include <vector>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SampleFilter.hpp"
#include "genfile/VariableInRangeSampleFilter.hpp"
#include "genfile/VariableInSetSampleFilter.hpp"

namespace genfile {
	SampleFilter::SampleFilter() {}
	SampleFilter::~SampleFilter() {}

	void SampleFilter::compute_failed_samples( genfile::CohortIndividualSource const& source, boost::function< void( std::size_t ) > callback ) const {
		std::vector< std::size_t > result ;
		for( std::size_t i = 0; i < source.get_number_of_individuals(); ++i ) {
			if( !test( source, i )) {
				callback( i ) ;
			}
		}
	}

	std::ostream& operator<< ( std::ostream& oStream, SampleFilter const& filter ) {
		filter.summarise( oStream ) ;
		return oStream ;
	}


	// Factory function for conditions.
	SampleFilter::UniquePtr SampleFilter::create( std::string const& spec ) {

		using namespace genfile::string_utils ;

		// for now just parse variable = value.
		std::vector< std::string > elts = split_and_strip( spec, "=", " \t\n" ) ;

		if( elts.size() !=2  ) {
			throw genfile::BadArgumentError( "condition_factory()", "spec=\"" + spec + "\"" ) ;
		}

		// remove quotes
		if( elts[1].size() >= 2 &&
			(
				( elts[1][0] == '"' && elts[1][elts[1].size()-1] == '"' )
				|| ( elts[1][0] == '\'' && elts[1][elts[1].size()-1] == '\'' )
			)
		) {
			elts[1] = elts[1].substr( 1, elts[1].size() - 2 ) ;
		}

		VariableInSetSampleFilter::UniquePtr filter( new VariableInSetSampleFilter( elts[0] ) ) ; 
		filter->add_level( elts[1] ) ;
		return SampleFilter::UniquePtr( filter.release() ) ;
	}
}
