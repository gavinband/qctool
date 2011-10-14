#include <string>
#include <boost/signals2/signal.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "AlleleFrequencyTestCallComparer.hpp"
#include "PairwiseCallComparer.hpp"

PairwiseCallComparer::UniquePtr PairwiseCallComparer::create( std::string const& spec ) {
	PairwiseCallComparer::UniquePtr result ;
	if( spec != "AlleleFrequencyTestCallComparer" ) {
		throw genfile::BadArgumentError( "PairwiseCallComparer::create()", "spec=\"" + spec + "\"" ) ;
	}
	result.reset( new AlleleFrequencyTestCallComparer() ) ;
	return result ;
}

