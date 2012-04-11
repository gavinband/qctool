
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include "genfile/Error.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "SampleBySampleComputation.hpp"
#include "RelatednessBayesFactorComputation.hpp"
#include "ConcordanceComputation.hpp"
#include "KinshipCoefficientComputation.hpp"

SampleBySampleComputation::UniquePtr SampleBySampleComputation::create( std::string const& name, appcontext::OptionProcessor const& options, appcontext::UIContext& ui_context ) {
	SampleBySampleComputation::UniquePtr result ;
	if( name == "relatedness" ) {
		result.reset(
			new RelatednessBayesFactorComputation(
				options,
				ui_context
			)
		) ;
	}
	else if( name == "kinship" ) {
		result.reset(
			new KinshipCoefficientComputation(
				options,
				ui_context
			)
		) ;
	}
	else if( name == "concordance" ) {
		result.reset(
			new ConcordanceComputation(
				options,
				ui_context
			)
		) ;
	}
	else if( name == "pairwise-non-missing-count" ) {
		result.reset(
			new PairwiseNonMissingnessComputation(
				options,
				ui_context
			)
		) ;
	}
	else {
		throw genfile::BadArgumentError( "SampleBySampleComputation::create()", "name=\"" + name + "\"" ) ;
	}
	return result ;
}

