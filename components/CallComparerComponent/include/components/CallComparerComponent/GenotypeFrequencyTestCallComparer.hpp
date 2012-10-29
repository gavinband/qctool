
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_Genotype_FREQUENCY_TEST_CALL_COMPARER_HPP
#define QCTOOL_Genotype_FREQUENCY_TEST_CALL_COMPARER_HPP

#include <string>
#include <map>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/function.hpp>
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "PairwiseCallComparer.hpp"

struct GenotypeFrequencyTestCallComparer: public PairwiseCallComparer {

	GenotypeFrequencyTestCallComparer() ;

	void compare(
		Eigen::MatrixXd const& left,
		Eigen::MatrixXd const& right,
		Callback
	) const ;
	
private:
	double m_threshhold ;
	boost::math::chi_squared_distribution< double > m_chi_squared ;	
} ;

#endif
