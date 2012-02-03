#ifndef QCTOOL_ALLELE_FREQUENCY_TEST_CALL_COMPARER_HPP
#define QCTOOL_ALLELE_FREQUENCY_TEST_CALL_COMPARER_HPP

#include <string>
#include <map>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/function.hpp>
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "PairwiseCallComparer.hpp"

struct AlleleFrequencyTestCallComparer: public PairwiseCallComparer {

	AlleleFrequencyTestCallComparer() ;

	std::map< std::string, genfile::VariantEntry > compare(
		genfile::SingleSNPGenotypeProbabilities const& left,
		genfile::SingleSNPGenotypeProbabilities const& right
	) const ;
	
private:
	double m_threshhold ;
	boost::math::chi_squared_distribution< double > m_chi_squared ;	
} ;

#endif
