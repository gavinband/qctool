
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <boost/signals2/signal.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "components/CallComparerComponent/AlleleFrequencyTestCallComparer.hpp"
#include "components/CallComparerComponent/PairwiseCallComparer.hpp"

PairwiseCallComparer::UniquePtr PairwiseCallComparer::create( std::string const& model ) {
	PairwiseCallComparer::UniquePtr result ;
	if( model == "AlleleFrequencyTest" ) {
		result.reset( new AlleleFrequencyTestCallComparer() ) ;
	}
	else if( model == "AcceptAll" ) {
		// empty.
	}
	else {
		throw genfile::BadArgumentError( "PairwiseCallComparer::create()", "model=\"" + model + "\"" ) ;
	}
	return result ;
}

