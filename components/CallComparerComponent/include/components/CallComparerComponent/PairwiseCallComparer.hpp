#ifndef QCTOOL_PAIRWISE_CALL_COMPARER_HPP
#define QCTOOL_PAIRWISE_CALL_COMPARER_HPP

#include <string>
#include <boost/signals2/signal.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/VariantEntry.hpp"

struct PairwiseCallComparer {
	typedef std::auto_ptr< PairwiseCallComparer > UniquePtr ;
	static UniquePtr create( std::string const& spec ) ;

	virtual ~PairwiseCallComparer() {}
	virtual std::map< std::string, genfile::VariantEntry > compare(
		genfile::SingleSNPGenotypeProbabilities const& left,
		genfile::SingleSNPGenotypeProbabilities const& right
	) const = 0 ;
} ;

#endif